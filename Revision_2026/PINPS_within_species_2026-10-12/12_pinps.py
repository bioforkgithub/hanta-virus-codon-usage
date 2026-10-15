#!/usr/bin/env python3
"""
Phase 2 step 3: within-species piN and piS (Nei & Gojobori 1986 counting) from the
codon alignments written by 11_align.py.

DEFINITIONS USED (all stated so they can be disagreed with)
-----------------------------------------------------------
Per PAIR of sequences (i, j):
  valid codons = positions where BOTH sequences carry a gap-free, ambiguity-free,
                 non-stop sense codon (pairwise deletion).  The mean number of
                 codons actually compared is reported.
  S_ij = sum over valid codons of (S(c_i) + S(c_j)) / 2 ;  N_ij likewise
      where S(c) = sum over the 3 codon positions of
            (number of synonymous single-nucleotide alternatives at that position)
          / (number of NON-STOP single-nucleotide alternatives at that position)
      and N(c) = 3 - S(c).  Stop codons are excluded from the site count (the
      MEGA / NG86 convention).  The variant that divides by 3 unconditionally is
      also computed and reported as piN_piS_stopexcl.
  Sd_ij, Nd_ij = synonymous / nonsynonymous differences, averaged over all
      shortest mutational pathways; pathways through a stop codon are discarded
      (if every pathway passes through one, all are used).
  pS_ij = Sd_ij / S_ij ,  pN_ij = Nd_ij / N_ij

Species-level estimates (two estimators, both reported):
  PRIMARY  piS = mean over all n(n-1)/2 pairs of pS_ij ; piN likewise.
           piN/piS = (mean pN) / (mean pS)  -- ratio of averages, not average of
           ratios (per-pair ratios are undefined whenever pS_ij = 0, which is
           common within a species).
  TOTALS   piS_tot = (sum over pairs Sd) / (sum over pairs S), piN_tot likewise.
           This estimator is additive over codons, which is what makes the codon
           bootstrap below exact; it is reported so the reader can see the two
           agree.
  Jukes-Cantor: dS_ij = -3/4 ln(1 - 4/3 pS_ij) applied ONLY where pS_ij < 0.75.
           piS_jc is the mean over the pairs where it is defined and the number of
           undefined pairs is reported.  max pS is reported per species so
           saturation can be checked directly -- the whole point of moving off
           between-species dN/dS.

CONFIDENCE INTERVALS
  "full" rows: CODON bootstrap (resample the L codon columns with replacement,
  1000 reps) on the additive TOTALS estimator.  A sequence bootstrap is NOT used:
  resampling sequences with replacement creates self-pairs with pS = pN = 0 and
  biases pi downward.  Sampling variance in the number of sequences is instead
  addressed by the equal-n subsampling rows, which is the control that matters
  here.
  subsample rows: 1000 draws of n sequences WITHOUT replacement; mean and
  2.5 / 97.5 percentiles of piN, piS and piN/piS across draws.

SEQUENCE SETS (row subsets of one alignment per species)
  lato_all   every usable complete N CDS under the [Organism] subtree
  ss_all     sensu stricto (separately named foreign lineages removed), all records
  ss_uniq    sensu stricto, one record per distinct nucleotide haplotype
  ss_thin    sensu stricto, one record per country x collection-year cell
  ss_animal  sensu stricto, records with a non-human /host
  clust95    largest single-linkage cluster at 95% nucleotide identity

Usage:  python3 -I 12_pinps.py <ANALYSIS_ROOT>
RNG: numpy PCG64, master seed 20261012, stream
     np.random.default_rng([20261012, species_index, set_index, n]) where n = 0
     for the codon bootstrap.  No wall-clock input anywhere.
"""
import sys, os, json, itertools
import numpy as np

MASTER_SEED = 20261012
N_REPS = 1000
N_BOOT = 1000
LADDER = [10, 12, 15, 20, 28, 47]

BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS, CODON_AA = [], {}
_i = 0
for b1 in BASES:
    for b2 in BASES:
        for b3 in BASES:
            CODONS.append(b1 + b2 + b3)
            CODON_AA[b1 + b2 + b3] = AAS[_i]
            _i += 1
SENSE = [c for c in CODONS if CODON_AA[c] != "*"]
CIDX = {c: k for k, c in enumerate(SENSE)}
NC = len(SENSE)
SENT = NC


def site_counts(exclude_stops=True):
    S = np.zeros(NC + 1)
    for c in SENSE:
        tot = 0.0
        for pos in range(3):
            alts = [c[:pos] + b + c[pos + 1:] for b in "ACGT" if b != c[pos]]
            if exclude_stops:
                alts = [a for a in alts if CODON_AA[a] != "*"]
            if not alts:
                continue
            tot += sum(1 for a in alts if CODON_AA[a] == CODON_AA[c]) / len(alts)
        S[CIDX[c]] = tot
    return S


def pathway_diffs(exclude_stop_pathways=False):
    SD = np.zeros((NC + 1, NC + 1))
    ND = np.zeros((NC + 1, NC + 1))
    for c1 in SENSE:
        for c2 in SENSE:
            diff = [p for p in range(3) if c1[p] != c2[p]]
            if not diff:
                continue
            if len(diff) == 1:
                syn = 1.0 if CODON_AA[c1] == CODON_AA[c2] else 0.0
                SD[CIDX[c1], CIDX[c2]] = syn
                ND[CIDX[c1], CIDX[c2]] = 1.0 - syn
                continue
            ok, allp = [], []
            for order in itertools.permutations(diff):
                cur, s, n, bad = c1, 0.0, 0.0, False
                for p in order:
                    nxt = cur[:p] + c2[p] + cur[p + 1:]
                    if CODON_AA[nxt] == "*":
                        bad = True
                    if CODON_AA[nxt] == CODON_AA[cur]:
                        s += 1
                    else:
                        n += 1
                    cur = nxt
                allp.append((s, n))
                if not bad:
                    ok.append((s, n))
            use = (ok if ok else allp) if exclude_stop_pathways else allp
            SD[CIDX[c1], CIDX[c2]] = sum(p[0] for p in use) / len(use)
            ND[CIDX[c1], CIDX[c2]] = sum(p[1] for p in use) / len(use)
    return SD, ND


def read_fasta(path):
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip().upper())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def encode(seqs):
    n = len(seqs)
    L = len(seqs[0]) // 3
    A = np.full((n, L), SENT, dtype=np.int16)
    NTM = np.full((n, 3 * L), 4, dtype=np.int8)
    b2i = {"A": 0, "C": 1, "G": 2, "T": 3}
    for i, s in enumerate(seqs):
        for j in range(L):
            k = CIDX.get(s[3 * j:3 * j + 3])
            if k is not None:
                A[i, j] = k
        for p, ch in enumerate(s):
            if ch in b2i:
                NTM[i, p] = b2i[ch]
    return A, NTM


def pairwise(A, NTM, SS, NS, SD, ND):
    """Full pairwise pN, pS, nucleotide p-distance, codons compared."""
    n, L = A.shape
    V = A != SENT
    Asafe = np.where(V, A, 0)
    Vn = NTM != 4
    pN = np.full((n, n), np.nan)
    pS = np.full((n, n), np.nan)
    ncmp = np.zeros((n, n))
    pdist = np.full((n, n), np.nan)
    for i in range(n - 1):
        m = V[i][None, :] & V[i + 1:]
        ai, bj = Asafe[i], Asafe[i + 1:]
        sd = (SD[ai[None, :], bj] * m).sum(1)
        nd = (ND[ai[None, :], bj] * m).sum(1)
        ss = ((SS[ai][None, :] + SS[bj]) * 0.5 * m).sum(1)
        ns = ((NS[ai][None, :] + NS[bj]) * 0.5 * m).sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            pS[i, i + 1:] = pS[i + 1:, i] = np.where(ss > 0, sd / ss, np.nan)
            pN[i, i + 1:] = pN[i + 1:, i] = np.where(ns > 0, nd / ns, np.nan)
        ncmp[i, i + 1:] = ncmp[i + 1:, i] = m.sum(1)
        mn = Vn[i][None, :] & Vn[i + 1:]
        tot = mn.sum(1)
        dif = ((NTM[i][None, :] != NTM[i + 1:]) & mn).sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            pdist[i, i + 1:] = pdist[i + 1:, i] = np.where(tot > 0, dif / tot, np.nan)
    return pN, pS, ncmp, pdist


def percodon_totals(A, idx, SS, NS, SD, ND):
    """Per-codon-column sums of Sd, S, Nd, N over all pairs inside idx."""
    sub = A[np.asarray(idx), :]
    n, L = sub.shape
    V = sub != SENT
    Asafe = np.where(V, sub, 0)
    T = np.zeros((4, L))
    for i in range(n - 1):
        m = V[i][None, :] & V[i + 1:]
        ai, bj = Asafe[i], Asafe[i + 1:]
        T[0] += (SD[ai[None, :], bj] * m).sum(0)
        T[1] += ((SS[ai][None, :] + SS[bj]) * 0.5 * m).sum(0)
        T[2] += (ND[ai[None, :], bj] * m).sum(0)
        T[3] += ((NS[ai][None, :] + NS[bj]) * 0.5 * m).sum(0)
    return T


def jc(p):
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p) & (p < 0.75)
    out[ok] = -0.75 * np.log(1.0 - (4.0 / 3.0) * p[ok])
    return out


def uptri(M, v):
    sub = M[np.ix_(v, v)]
    iu = np.triu_indices(len(v), 1)
    return sub[iu]


def seg_sites(NTM, A, idx):
    sub = NTM[np.asarray(idx), :]
    nseg = 0
    for p in range(sub.shape[1]):
        col = sub[:, p]
        col = col[col != 4]
        if len(col) > 1 and len(np.unique(col)) > 1:
            nseg += 1
    Asub = A[np.asarray(idx), :]
    poly = aach = 0
    for j in range(Asub.shape[1]):
        col = Asub[:, j]
        col = col[col != SENT]
        u = np.unique(col)
        if len(u) > 1:
            poly += 1
            if len({CODON_AA[SENSE[k]] for k in u}) > 1:
                aach += 1
    return nseg, poly, aach


def clusters_95(pdist, thresh=0.05):
    n = pdist.shape[0]
    adj = np.isfinite(pdist) & (pdist <= thresh)
    np.fill_diagonal(adj, True)
    seen = np.zeros(n, bool)
    comps = []
    for s in range(n):
        if seen[s]:
            continue
        stack, comp = [s], []
        seen[s] = True
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in np.where(adj[u] & ~seen)[0]:
                seen[v] = True
                stack.append(v)
        comps.append(sorted(comp))
    comps.sort(key=lambda c: (-len(c), c[0]))
    return comps


def prune_to_diameter(pdist, members, thresh=0.05):
    """Greedy prune of a single-linkage cluster down to DIAMETER <= thresh.

    Single linkage chains: two sequences 9% apart can sit in the same 95%
    cluster via intermediates.  This removes, one at a time, the member
    involved in the most violating pairs (ties broken by the larger sum of
    distances, then by index, so the result is deterministic) until no pair in
    the set exceeds thresh.  The result is a set whose every pairwise distance
    is within the threshold, i.e. a defensible single-population sample.
    """
    cur = list(members)
    while len(cur) > 2:
        sub = pdist[np.ix_(cur, cur)]
        bad = np.isfinite(sub) & (sub > thresh)
        np.fill_diagonal(bad, False)
        if not bad.any():
            break
        cnt = bad.sum(1)
        tot = np.where(np.isfinite(sub), sub, 0.0).sum(1)
        worst = max(range(len(cur)), key=lambda k: (cnt[k], tot[k], -cur[k]))
        cur.pop(worst)
    return sorted(cur)


def main():
    root = sys.argv[1]
    aldir = os.path.join(root, "results", "alignments")
    resdir = os.path.join(root, "results")

    # PRIMARY: Nei & Gojobori (1986) as originally specified and as implemented
    # in Biopython 1.83 cal_dn_ds(method="NG86") -- site counts divide by 3 with
    # stop-codon neighbours counted as nonsynonymous, and codon-pair differences
    # average over ALL shortest mutational pathways.  Verified to agree with
    # Biopython to 1.8e-15 relative on 25 random 200-codon pairs
    # (scripts/13_validate_ng86.py).
    SS = site_counts(False)
    NS = np.where(np.arange(NC + 1) < NC, 3.0 - SS, 0.0)
    SD, ND = pathway_diffs(exclude_stop_pathways=False)
    # ROBUSTNESS variant (MEGA convention): stop-codon changes excluded from the
    # site denominators and stop-passing pathways discarded.
    SS0 = site_counts(True)
    NS0 = np.where(np.arange(NC + 1) < NC, 3.0 - SS0, 0.0)
    SD0, ND0 = pathway_diffs(exclude_stop_pathways=True)

    meta = {}
    with open(os.path.join(resdir, "sets", "sequence_metadata.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            meta[(r["abbrev"], r["accession"])] = r

    order = ["HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "CHOV",
             "TULV", "NVAV", "MJNV", "TPMV", "SWSV", "CBNV", "JJUV",
             "RKPV", "PHV", "ASAV"]
    setnames = ["lato_all", "ss_all", "ss_uniq", "ss_thin", "ss_animal",
                "clust95", "clust95d"]

    rows, diag = [], {}
    for si, ab in enumerate(order):
        path = os.path.join(aldir, "%s.codon.aln.fasta" % ab)
        if not os.path.exists(path):
            continue
        recs = read_fasta(path)
        names = [r[0] for r in recs]
        A, NTM = encode([r[1] for r in recs])
        pN, pS, ncmp, pdist = pairwise(A, NTM, SS, NS, SD, ND)
        cat = meta[(ab, names[0])]["category"]
        hc = meta[(ab, names[0])]["host_class"]

        iall = list(range(len(names)))
        iss = [i for i, nm in enumerate(names)
               if meta[(ab, nm)]["in_sensu_stricto"] == "True"]
        seen, iuniq = set(), []
        for i in iss:
            if recs[i][1] not in seen:
                seen.add(recs[i][1]); iuniq.append(i)
        cells, ithin = set(), []
        for i in iss:
            m = meta[(ab, names[i])]
            if m["country"] and m["year"]:
                key = (m["country"], m["year"])
                if key in cells:
                    continue
                cells.add(key)
            ithin.append(i)
        ianim = [i for i, nm in enumerate(names)
                 if meta[(ab, nm)]["in_sensu_stricto"] == "True"
                 and meta[(ab, nm)]["is_human"] != "True"
                 and meta[(ab, nm)]["host"].strip()]
        comps = clusters_95(pdist)
        iclust = comps[0] if comps else []
        # diameter-constrained single-population set: prune every single-linkage
        # cluster to diameter <= 5% and keep the largest survivor
        pruned = [prune_to_diameter(pdist, c) for c in comps if len(c) >= 2]
        pruned.sort(key=lambda c: (-len(c), c[0]))
        iclustd = pruned[0] if pruned else []

        idxmap = dict(lato_all=iall, ss_all=iss, ss_uniq=iuniq,
                      ss_thin=ithin, ss_animal=ianim, clust95=iclust,
                      clust95d=iclustd)
        diag[ab] = dict(
            category=cat, host_class=hc, aln_codons=int(A.shape[1]),
            n_lato=len(iall), n_ss=len(iss), n_uniq=len(iuniq),
            n_thin=len(ithin), n_animal=len(ianim),
            n_clusters_95=len(comps), cluster_sizes_95=[len(c) for c in comps[:12]],
            largest_cluster_internal_max_pdist=(
                float(np.nanmax(uptri(pdist, np.asarray(iclust))))
                if len(iclust) > 1 else None),
            clust95_accessions=[names[i] for i in iclust][:60],
            n_clust95d=len(iclustd),
            clust95d_diameter=(float(np.nanmax(uptri(pdist, np.asarray(iclustd))))
                               if len(iclustd) > 1 else None),
            clust95d_accessions=[names[i] for i in iclustd][:60],
            clust95d_countries=sorted({meta[(ab, names[i])]["country"]
                                       for i in iclustd
                                       if meta[(ab, names[i])]["country"]}),
            clust95d_hosts=sorted({meta[(ab, names[i])]["host"] for i in iclustd
                                   if meta[(ab, names[i])]["host"].strip()}))

        for sj, sn in enumerate(setnames):
            idx = np.asarray(idxmap[sn])
            if len(idx) < 2:
                continue
            an, as_ = uptri(pN, idx), uptri(pS, idx)
            fin = np.isfinite(an) & np.isfinite(as_)
            an, as_ = an[fin], as_[fin]
            ad = uptri(pdist, idx)
            ad = ad[np.isfinite(ad)]
            piN, piS_ = float(np.mean(an)), float(np.mean(as_))
            jn, js = jc(an), jc(as_)
            jn_ok, js_ok = jn[np.isfinite(jn)], js[np.isfinite(js)]
            nseg, poly, aach = seg_sites(NTM, A, idx)

            T = percodon_totals(A, idx, SS, NS, SD, ND)
            T0 = percodon_totals(A, idx, SS0, NS0, SD0, ND0)
            piS_tot = T[0].sum() / T[1].sum() if T[1].sum() > 0 else np.nan
            piN_tot = T[2].sum() / T[3].sum() if T[3].sum() > 0 else np.nan
            r0 = ((T0[2].sum() / T0[3].sum()) / (T0[0].sum() / T0[1].sum())
                  if T0[0].sum() > 0 else np.nan)

            L = A.shape[1]
            rng = np.random.default_rng([MASTER_SEED, si, sj, 0])
            bn, bs, br = np.empty(N_BOOT), np.empty(N_BOOT), np.empty(N_BOOT)
            for b in range(N_BOOT):
                col = rng.integers(0, L, L)
                sd, s_, nd, n_ = T[0][col].sum(), T[1][col].sum(), T[2][col].sum(), T[3][col].sum()
                bs[b] = sd / s_ if s_ > 0 else np.nan
                bn[b] = nd / n_ if n_ > 0 else np.nan
                br[b] = bn[b] / bs[b] if bs[b] > 0 else np.nan

            mset = [meta[(ab, names[i])] for i in idx]
            ncountry = len({m["country"] for m in mset if m["country"]})
            nhuman = sum(1 for m in mset if m["is_human"] == "True")
            piS_jc_v = float(np.mean(js_ok)) if len(js_ok) else np.nan
            sat = bool((piS_jc_v >= 0.75) or (float(np.max(as_)) >= 0.5)
                       or int(np.sum(as_ >= 0.75)) > 0)
            lineage_pool = bool((float(np.median(ad)) > 0.05 if len(ad) else False)
                                or (float(np.max(ad)) > 0.10 if len(ad) else False))
            rows.append(dict(
                abbrev=ab, category=cat, host_class=hc,
                category_tulv_moved=("pathogenic" if ab == "TULV" else cat),
                seqset=sn, n_level="full",
                n_sequences=len(idx), n_countries=ncountry, n_human=nhuman,
                n_reps=N_BOOT, resample="codon bootstrap",
                aln_codons=L,
                flag_piS_saturated=sat, flag_lineage_pool=lineage_pool,
                usable_for_group_test=bool(len(idx) >= 10 and not sat
                                           and not lineage_pool),
                mean_codons_compared=float(np.mean(uptri(ncmp, idx))),
                n_segregating_nt=nseg, n_polymorphic_codons=poly,
                n_aa_changing_codons=aach,
                piN=piN, piS=piS_, piN_piS=(piN / piS_ if piS_ > 0 else np.nan),
                piN_tot=float(piN_tot), piS_tot=float(piS_tot),
                piN_piS_tot=(float(piN_tot / piS_tot) if piS_tot > 0 else np.nan),
                piN_jc=float(np.mean(jn_ok)) if len(jn_ok) else np.nan,
                piS_jc=float(np.mean(js_ok)) if len(js_ok) else np.nan,
                piN_piS_jc=(float(np.mean(jn_ok) / np.mean(js_ok))
                            if len(js_ok) and np.mean(js_ok) > 0 else np.nan),
                piN_lo=float(np.nanpercentile(bn, 2.5)),
                piN_hi=float(np.nanpercentile(bn, 97.5)),
                piS_lo=float(np.nanpercentile(bs, 2.5)),
                piS_hi=float(np.nanpercentile(bs, 97.5)),
                ratio_mean=float(np.nanmean(br)),
                ratio_lo=float(np.nanpercentile(br, 2.5)),
                ratio_hi=float(np.nanpercentile(br, 97.5)),
                max_pS=float(np.max(as_)), max_pN=float(np.max(an)),
                median_pdist_nt=float(np.median(ad)) if len(ad) else np.nan,
                max_pdist_nt=float(np.max(ad)) if len(ad) else np.nan,
                n_pairs=int(len(an)),
                n_pairs_pS_zero=int(np.sum(as_ == 0)),
                n_pairs_pS_ge_075=int(np.sum(as_ >= 0.75)),
                n_pairs_jc_undefined_S=int(np.sum(~np.isfinite(js))),
                piN_piS_stopexcl=float(r0)))

            for n in LADDER:
                if len(idx) <= n:
                    continue
                rng = np.random.default_rng([MASTER_SEED, si, sj, n])
                vn = np.empty(N_REPS); vs = np.empty(N_REPS); vr = np.empty(N_REPS)
                for b in range(N_REPS):
                    v = idx[rng.choice(len(idx), n, replace=False)]
                    a_n, a_s = uptri(pN, v), uptri(pS, v)
                    f = np.isfinite(a_n) & np.isfinite(a_s)
                    mn = float(np.mean(a_n[f])); ms = float(np.mean(a_s[f]))
                    vn[b] = mn; vs[b] = ms
                    vr[b] = mn / ms if ms > 0 else np.nan
                rows.append(dict(
                    abbrev=ab, category=cat, host_class=hc,
                    category_tulv_moved=("pathogenic" if ab == "TULV" else cat),
                    seqset=sn, n_level=str(n),
                    n_sequences=n, n_reps=N_REPS,
                    flag_piS_saturated=sat, flag_lineage_pool=lineage_pool,
                    usable_for_group_test=bool(not sat and not lineage_pool),
                    resample="subsample without replacement", aln_codons=L,
                    piN=float(np.mean(vn)), piS=float(np.mean(vs)),
                    piN_piS=float(np.nanmean(vr)),
                    piN_lo=float(np.percentile(vn, 2.5)),
                    piN_hi=float(np.percentile(vn, 97.5)),
                    piS_lo=float(np.percentile(vs, 2.5)),
                    piS_hi=float(np.percentile(vs, 97.5)),
                    ratio_mean=float(np.nanmean(vr)),
                    ratio_lo=float(np.nanpercentile(vr, 2.5)),
                    ratio_hi=float(np.nanpercentile(vr, 97.5)),
                    n_pairs=n * (n - 1) // 2,
                    n_subsamples_ratio_undefined=int(np.sum(~np.isfinite(vr)))))
        print("%-5s n_lato=%-4d n_ss=%-4d n_uniq=%-4d n_thin=%-4d n_anim=%-4d "
              "clusters95=%-3d largest=%-4d diam-pruned=%d"
              % (ab, len(iall), len(iss), len(iuniq), len(ithin), len(ianim),
                 len(comps), len(iclust), len(iclustd)), flush=True)

    COLS = ["abbrev", "category", "category_tulv_moved", "host_class", "seqset",
            "n_level", "n_sequences", "n_countries", "n_human",
            "flag_piS_saturated", "flag_lineage_pool", "usable_for_group_test",
            "n_reps", "resample", "aln_codons", "mean_codons_compared",
            "n_segregating_nt", "n_polymorphic_codons", "n_aa_changing_codons",
            "piN", "piS", "piN_piS", "piN_tot", "piS_tot", "piN_piS_tot",
            "piN_jc", "piS_jc", "piN_piS_jc",
            "piN_lo", "piN_hi", "piS_lo", "piS_hi",
            "ratio_mean", "ratio_lo", "ratio_hi",
            "max_pS", "max_pN", "median_pdist_nt", "max_pdist_nt",
            "n_pairs", "n_pairs_pS_zero", "n_pairs_pS_ge_075",
            "n_pairs_jc_undefined_S", "n_subsamples_ratio_undefined",
            "piN_piS_stopexcl"]

    def fmt(v):
        if isinstance(v, float):
            return "NA" if not np.isfinite(v) else "%.6g" % v
        return str(v)

    with open(os.path.join(resdir, "pinps_per_species.tsv"), "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c, "")) for c in COLS) + "\n")
    with open(os.path.join(resdir, "pinps_diagnostics.json"), "w") as fh:
        json.dump(dict(species=diag, master_seed=MASTER_SEED, n_reps=N_REPS,
                       n_boot=N_BOOT, ladder=LADDER,
                       site_count_check=dict(
                           primary_convention="NG86 original / Biopython NG86:"
                           " sites /3, stop neighbours nonsynonymous, all"
                           " shortest pathways",
                           S_TTT_primary=float(SS[CIDX["TTT"]]),
                           S_TGG_primary=float(SS[CIDX["TGG"]]),
                           S_CGA_primary=float(SS[CIDX["CGA"]]),
                           S_CGA_stopexcl=float(SS0[CIDX["CGA"]]),
                           mean_S_over_61_sense_codons_primary=float(SS[:NC].mean()),
                           mean_S_over_61_sense_codons_stopexcl=float(SS0[:NC].mean()))),
                  fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
