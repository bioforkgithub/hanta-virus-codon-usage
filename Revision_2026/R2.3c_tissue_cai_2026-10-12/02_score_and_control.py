#!/usr/bin/env python3
"""R2.3c + R4.1 step 2 - score viral and host CDS against all 51 tissue
references, with a positive control that tests whether the method has any
tissue resolution at all.

Three questions, each answered with its own output:

A. R2.3c: is viral CAI higher under the lung and kidney references
   than under control tissues? -> rank of lung/kidney among 51 tissues, with a
   bootstrap CI over the reference gene set.

B. Validity (our own additional check): do genes that are themselves
   tissue-ENRICHED score highest under their own tissue's reference? If a
   lung-enriched gene set does not rank lung top, then no tissue-level CAI
   signal is recoverable for anything, virus included, and both R2.3c and R4.1
   are answered by that fact.

C. How different are the 51 references from each other at all?
"""
import collections
import csv
import glob
import json
import math
import os
import random
import re
import sys

from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cai_lib  # noqa: E402

HPA = sys.argv[1]
CDSFA = sys.argv[2]
W1 = sys.argv[3]
RESULTS = sys.argv[4]
LABELS = sys.argv[5]
TOPN = 500
BOOT = int(sys.argv[6]) if len(sys.argv) > 6 else 200
random.seed(20261012)

# ---------------------------------------------------------------- references
gene_re = re.compile(r"gene:(ENSG\d+)")
bio_re = re.compile(r"transcript_biotype:(\S+)")
sym_re = re.compile(r"gene_symbol:(\S+)")
best = {}


def consider(hdr, seq):
    g, b = gene_re.search(hdr), bio_re.search(hdr)
    if not g or not b or b.group(1) != "protein_coding":
        return
    seq = seq.upper()
    if len(seq) < 300 or len(seq) % 3 or set(seq) - set("ACGT"):
        return
    if not seq.startswith("ATG") or seq[-3:] not in cai_lib.STOPS:
        return
    cur = best.get(g.group(1))
    if cur is None or len(seq) > len(cur[1]):
        s = sym_re.search(hdr)
        best[g.group(1)] = (s.group(1) if s else "", seq)


hdr, buf = None, []
with open(CDSFA) as fh:
    for line in fh:
        if line.startswith(">"):
            if hdr:
                consider(hdr, "".join(buf))
            hdr, buf = line, []
        else:
            buf.append(line.strip())
if hdr:
    consider(hdr, "".join(buf))
gcounts = {g: cai_lib.codon_counts(s) for g, (_s, s) in best.items()}

expr = collections.defaultdict(dict)
with open(HPA) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    i_g, i_t, i_v = head.index("Gene"), head.index("Tissue"), head.index("nTPM")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) > max(i_g, i_t, i_v):
            try:
                expr[p[i_t]][p[i_g]] = float(p[i_v])
            except ValueError:
                pass
tissues = sorted(expr)

topgenes = {}
for t in tissues:
    r = sorted(((v, g) for g, v in expr[t].items() if g in gcounts), reverse=True)
    topgenes[t] = [g for _v, g in r[:TOPN]]


def weights_from_genes(genes):
    tot = collections.Counter()
    for g in genes:
        tot.update(gcounts[g])
    return cai_lib.weights_proper(cai_lib.fractions_from_counts(tot))


W = {t: weights_from_genes(topgenes[t]) for t in tissues}

# ------------------------------------------------------------- sequence sets
labels = {}
with open(LABELS) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        labels[row["accession"]] = row

N_PAT = re.compile(r"nucleocapsid|nucleoprotein", re.I)
GPC_PAT = re.compile(r"glycoprotein|precursor|\bGPC\b|Gn|Gc", re.I)
L_PAT = re.compile(r"polymerase|\bL protein\b", re.I)

sets = collections.defaultdict(dict)


def clean(seq):
    s = str(seq).upper()
    if len(s) % 3 or set(s) - set("ACGT") or len(s) < 300:
        return None
    return s


for p in sorted(glob.glob(os.path.join(W1, "Genus_Diversity", "gbk", "*.gbk"))):
    for rec in SeqIO.parse(p, "genbank"):
        if rec.id not in labels:
            continue  # keep exactly the manuscript's 22 species
        ab = labels[rec.id]["abbrev"]
        for f in rec.features:
            if f.type != "CDS":
                continue
            prod = (f.qualifiers.get("product") or [""])[0]
            if N_PAT.search(prod):
                s = clean(f.extract(rec.seq))
                if s:
                    sets["genus22_N"]["%s_%s" % (ab, rec.id)] = s

for p in sorted(glob.glob(os.path.join(W1, "Epuyen_outbreak", "gbk", "*"))):
    for rec in SeqIO.parse(p, "genbank"):
        for f in rec.features:
            if f.type not in ("CDS", "RNA-dependent RNA polymerase"):
                continue
            prod = (f.qualifiers.get("product") or [""])[0] or f.type
            s = clean(f.extract(rec.seq))
            if not s:
                continue
            if N_PAT.search(prod):
                sets["epuyen_ANDV_N"][rec.id] = s
            elif L_PAT.search(prod):
                sets["epuyen_ANDV_L"][rec.id] = s
            elif GPC_PAT.search(prod):
                sets["epuyen_ANDV_GPC"][rec.id] = s

for nm, path in (("lung_HKG_panel", os.path.join(W1, "Lung_HKG", "lung_hkg_cds.fasta")),
                 ("kidney_HKG_panel", os.path.join(W1, "Kidney_HKG", "kidney_hkg_cds.fasta"))):
    if os.path.exists(path):
        for rec in SeqIO.parse(path, "fasta"):
            s = clean(rec.seq)
            if s:
                sets[nm][rec.id] = s

# positive control: tissue-ENRICHED gene sets (max in tissue, >=5x mean of rest)
enriched = {}
for t in tissues:
    out = []
    for g in gcounts:
        v = expr[t].get(g, 0.0)
        if v < 10:
            continue
        others = [expr[o].get(g, 0.0) for o in tissues if o != t]
        m = sum(others) / len(others)
        if v >= 5 * max(m, 1e-9) and v == max([v] + others):
            out.append((v, g))
    out.sort(reverse=True)
    enriched[t] = [g for _v, g in out[:100]]
for t in ("lung", "kidney", "liver", "testis", "skeletal muscle", "pancreas"):
    if enriched.get(t):
        sets["ENRICHED_" + re.sub(r"\W+", "_", t)] = {
            (best[g][0] or g): best[g][1] for g in enriched[t]}

# ------------------------------------------------------------------- scoring
def score_set(seqs, w):
    vals = []
    for s in seqs.values():
        c, used, _ = cai_lib.cai_proper(s, w)
        if c is not None:
            vals.append(c)
    return vals


def mean(v):
    return sum(v) / len(v)


summary = {}
for name, seqs in sorted(sets.items()):
    per_tissue = {t: mean(score_set(seqs, W[t])) for t in tissues}
    order = sorted(per_tissue, key=lambda t: -per_tissue[t])
    rank = {t: i + 1 for i, t in enumerate(order)}
    with open(os.path.join(RESULTS, "cai_by_tissue_%s.tsv" % name), "w") as fh:
        fh.write("tissue\tmean_cai\trank\n")
        for t in order:
            fh.write("%s\t%.6f\t%d\n" % (t, per_tissue[t], rank[t]))
    vals = sorted(per_tissue.values())
    summary[name] = {
        "n_seqs": len(seqs),
        "lung_rank": rank.get("lung"), "kidney_rank": rank.get("kidney"),
        "n_tissues": len(tissues),
        "lung_cai": per_tissue.get("lung"), "kidney_cai": per_tissue.get("kidney"),
        "best_tissue": order[0], "best_cai": per_tissue[order[0]],
        "worst_tissue": order[-1], "worst_cai": per_tissue[order[-1]],
        "median_cai": vals[len(vals) // 2],
        "spread_max_minus_min": vals[-1] - vals[0],
        "top5": order[:5], "bottom5": order[-5:],
    }

# per-sequence matrix for the two sets reported
for name in ("genus22_N", "epuyen_ANDV_N"):
    if name not in sets:
        continue
    with open(os.path.join(RESULTS, "cai_matrix_%s.tsv" % name), "w") as fh:
        fh.write("sequence\t" + "\t".join(tissues) + "\n")
        for sid, s in sorted(sets[name].items()):
            row = []
            for t in tissues:
                c, _u, _k = cai_lib.cai_proper(s, W[t])
                row.append("%.6f" % c if c is not None else "NA")
            fh.write(sid + "\t" + "\t".join(row) + "\n")

# -------------------------------------------- bootstrap over reference genes
boot_out = {}
for name in ("genus22_N", "epuyen_ANDV_N"):
    if name not in sets:
        continue
    seqs = sets[name]
    ranks_lung, ranks_kidney = [], []
    for _ in range(BOOT):
        pt = {}
        for t in tissues:
            gs = [random.choice(topgenes[t]) for _ in topgenes[t]]
            pt[t] = mean(score_set(seqs, weights_from_genes(gs)))
        order = sorted(pt, key=lambda t: -pt[t])
        ranks_lung.append(order.index("lung") + 1)
        ranks_kidney.append(order.index("kidney") + 1)
    ranks_lung.sort()
    ranks_kidney.sort()
    boot_out[name] = {
        "reps": BOOT,
        "lung_rank_median": ranks_lung[BOOT // 2],
        "lung_rank_ci95": [ranks_lung[int(0.025 * BOOT)], ranks_lung[int(0.975 * BOOT)]],
        "kidney_rank_median": ranks_kidney[BOOT // 2],
        "kidney_rank_ci95": [ranks_kidney[int(0.025 * BOOT)], ranks_kidney[int(0.975 * BOOT)]],
    }

# --------------------------------------- how distinguishable are references?
codons = sorted(set().union(*[set(w) for w in W.values()]))
spread = []
for c in codons:
    vs = [W[t].get(c) for t in tissues if c in W[t]]
    spread.append((max(vs) - min(vs), c, min(vs), max(vs)))
spread.sort(reverse=True)
dists = []
for i, a in enumerate(tissues):
    for b in tissues[i + 1:]:
        shared = [c for c in codons if c in W[a] and c in W[b]]
        d = math.sqrt(sum((W[a][c] - W[b][c]) ** 2 for c in shared) / len(shared))
        dists.append((d, a, b))
dists.sort()
ref_sep = {
    "n_tissues": len(tissues),
    "codons_most_variable_across_tissues": [
        {"codon": c, "aa": cai_lib.CODON2AA.get(c), "min_w": round(lo, 4),
         "max_w": round(hi, 4), "range": round(r, 4)} for r, c, lo, hi in spread[:10]],
    "codons_least_variable": [
        {"codon": c, "range": round(r, 4)} for r, c, lo, hi in spread[-5:]],
    "rms_weight_distance_between_tissue_pairs": {
        "min": dists[0][0], "min_pair": [dists[0][1], dists[0][2]],
        "median": dists[len(dists) // 2][0],
        "max": dists[-1][0], "max_pair": [dists[-1][1], dists[-1][2]],
    },
    "lung_vs_kidney_rms_distance": next(
        d for d, a, b in dists if {a, b} == {"lung", "kidney"}),
}

out = {"topN_reference_genes": TOPN, "sets": summary,
       "bootstrap_tissue_rank": boot_out, "reference_separation": ref_sep}
with open(os.path.join(RESULTS, "R2.3c_summary.json"), "w") as fh:
    json.dump(out, fh, indent=2, sort_keys=True)

print("%-26s %5s %9s %9s %10s %10s %-20s %-20s" % (
    "set", "n", "lungRank", "kidRank", "lungCAI", "kidCAI", "best tissue", "worst tissue"))
for k, v in sorted(summary.items()):
    print("%-26s %5d %9s %9s %10.5f %10.5f %-20s %-20s" % (
        k, v["n_seqs"], v["lung_rank"], v["kidney_rank"], v["lung_cai"],
        v["kidney_cai"], v["best_tissue"][:20], v["worst_tissue"][:20]))
print()
print(json.dumps(boot_out, indent=2))
print()
print(json.dumps(ref_sep["rms_weight_distance_between_tissue_pairs"], indent=2))
print("lung vs kidney RMS weight distance: %.5f" % ref_sep["lung_vs_kidney_rms_distance"])
