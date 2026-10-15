#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""C3 decisive check: re-derive NG86 ps/pn from the raw codon alignment with a
from-scratch implementation (own codon table, own synonymous-site counter, own
equal-weight pathway counter), and ask three questions:

  Q1 does it reproduce the manuscript's ps and pn exactly?
  Q2 is the manuscript's dNdS column pn/ps, or is it the Jukes-Cantor-corrected
     dN/dS that the paper's Methods claim ("Nei & Gojobori 1986")?
  Q3 for how many pairs would the real NG86 dN/dS be computable at all?

Nothing is imported from Work1. Work1 is read-only here.
"""
import csv, itertools, math
from pathlib import Path
from statistics import mean

GD = Path("/PATH/TO/"
          "Project_Hantavirus/Work1/Genus_Diversity")
A = Path("/PATH/TO/"
         "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12")

B = "TCAG"
AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
TBL = {}
k = 0
for a in B:
    for b in B:
        for c in B:
            pass
# build the standard table in the canonical TCAG x TCAG x TCAG order
TBL = {}
i = 0
for a in B:
    for b in B:
        for c in B:
            TBL[a + b + c] = AAS[i]; i += 1
assert TBL["ATG"] == "M" and TBL["TGA"] == "*" and TBL["TTT"] == "F"
assert TBL["AAA"] == "K" and TBL["GGG"] == "G" and TBL["TGG"] == "W"


def syn_sites(cod):
    """NG86 synonymous sites in one codon: sum over 3 positions of (#synonymous
    changes at that position)/3."""
    s = 0.0
    aa = TBL[cod]
    for p in range(3):
        for nt in "TCAG":
            if nt == cod[p]:
                continue
            alt = cod[:p] + nt + cod[p + 1:]
            if TBL[alt] == aa:
                s += 1.0 / 3.0
    return s


def path_diff(c1, c2):
    """Equal-weighting over all shortest mutational pathways; returns (Sd, Nd)."""
    diff = [p for p in range(3) if c1[p] != c2[p]]
    if not diff:
        return 0.0, 0.0
    sd = nd = 0.0
    npaths = 0
    for perm in itertools.permutations(diff):
        cur = c1
        s = n = 0
        ok = True
        for p in perm:
            nxt = cur[:p] + c2[p] + cur[p + 1:]
            if TBL[nxt] == "*":          # pathways through a stop are excluded
                ok = False
                break
            if TBL[nxt] == TBL[cur]:
                s += 1
            else:
                n += 1
            cur = nxt
        if ok:
            sd += s; nd += n; npaths += 1
    if npaths == 0:
        return 0.0, float(len(diff))
    return sd / npaths, nd / npaths


def jc(p):
    if p >= 0.75:
        return None
    return -0.75 * math.log(1 - 4.0 * p / 3.0)


# ---- read alignment
seqs = {}
name = None
for ln in (GD / "genus_diversity_codon_aligned.fasta").read_text().splitlines():
    if ln.startswith(">"):
        name = ln[1:].split("|")[0]
        seqs[name] = []
    elif ln.strip():
        seqs[name].append(ln.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
L = len(next(iter(seqs.values())))
print(f"alignment: {len(seqs)} seqs, {L} nt = {L//3} codons")

codons = {k: [v[i:i + 3] for i in range(0, L, 3)] for k, v in seqs.items()}

ms = {}
for r in csv.DictReader((GD / "results" / "pairwise_dnds_matrix.tsv").open(), delimiter="\t"):
    ms[frozenset((r["species_A"], r["species_B"]))] = r

labels = [r["abbrev"] for r in csv.DictReader((A / "inputs" / "labels.tsv").open(), delimiter="\t")]
cat = {r["abbrev"]: r["category"] for r in
       csv.DictReader((A / "inputs" / "labels.tsv").open(), delimiter="\t")}

worst_ps = worst_pn = worst_ncod = 0.0
worst_dnds_vs_pnps = 0.0
worst_dnds_vs_realNG = 0.0
real_defined = 0
rows = []
for a, b in itertools.combinations(labels, 2):
    ca, cb = codons[a], codons[b]
    S = N = Sd = Nd = 0.0
    nc = 0
    for c1, c2 in zip(ca, cb):
        if "-" in c1 or "-" in c2:
            continue
        nc += 1
        if TBL.get(c1, "*") == "*" or TBL.get(c2, "*") == "*":
            continue
        s1, s2 = syn_sites(c1), syn_sites(c2)
        S += (s1 + s2) / 2
        N += 3 - (s1 + s2) / 2
        sd, nd = path_diff(c1, c2)
        Sd += sd; Nd += nd
    ps, pn = Sd / S, Nd / N
    r = ms[frozenset((a, b))]
    worst_ps = max(worst_ps, abs(ps - float(r["ps"])))
    worst_pn = max(worst_pn, abs(pn - float(r["pn"])))
    worst_ncod = max(worst_ncod, abs(nc - int(r["n_codons_compared"])))
    # Q2: is the published column pn/ps, or the JC-corrected ratio?
    worst_dnds_vs_pnps = max(worst_dnds_vs_pnps, abs(float(r["dNdS"]) - pn / ps))
    dS, dN = jc(ps), jc(pn)
    if dS is not None and dN is not None and dS > 0:
        real_defined += 1
        worst_dnds_vs_realNG = max(worst_dnds_vs_realNG, abs(float(r["dNdS"]) - dN / dS))
    rows.append((a, b, cat[a], cat[b], ps, pn, pn / ps,
                 dN, dS, (dN / dS) if (dS and dN) else None))

print(f"\nQ1  from-scratch NG86 vs the manuscript's own table (231 pairs):")
print(f"    max |ps_recomputed  - ps_published |  = {worst_ps:.3e}")
print(f"    max |pn_recomputed  - pn_published |  = {worst_pn:.3e}")
print(f"    max |n_codons diff|                   = {worst_ncod:.0f}")
print(f"\nQ2  what is the published dNdS column?")
print(f"    max |dNdS_published - pn/ps        |  = {worst_dnds_vs_pnps:.3e}   <- pure rounding")
print(f"    max |dNdS_published - JC dN/dS     |  = {worst_dnds_vs_realNG:.3e}   "
      f"(over the {real_defined} pairs where JC dN/dS exists at all)")
print(f"\nQ3  pairs where a genuine NG86 (JC-corrected) dN/dS is defined: "
      f"{real_defined}/231")

# what the real NG86 group means would be, if anyone asked
GROUPS = ["pathogenic-pathogenic", "cross", "environmental-environmental"]
def g_of(ca, cb):
    if ca == cb == "pathogenic": return GROUPS[0]
    if ca == cb == "environmental": return GROUPS[2]
    return "cross"
print("\n    group means of the from-scratch JC-corrected dN/dS, where defined:")
for g in GROUPS:
    v = [r[9] for r in rows if g_of(r[2], r[3]) == g and r[9] is not None]
    tot = sum(1 for r in rows if g_of(r[2], r[3]) == g)
    print(f"      {g:30s} n={len(v):3d}/{tot:3d} mean={mean(v):.4f}" if v
          else f"      {g:30s} n=0/{tot}")
