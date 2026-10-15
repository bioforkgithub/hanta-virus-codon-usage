#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""C8: are the three estimators really on the SAME data?  C9: is the NG86
'ordering preserved' result driven by non-random dropout?"""
import csv, itertools, math
from pathlib import Path
from statistics import mean, median
from scipy import stats as sps

GD = Path("/PATH/TO/"
          "Project_Hantavirus/Work1/Genus_Diversity")
A = Path("/PATH/TO/"
         "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12")

seqs, nm = {}, None
for ln in (A / "inputs" / "S22_codon_aln.fasta").read_text().splitlines():
    if ln.startswith(">"):
        nm = ln[1:].strip(); seqs[nm] = []
    elif ln.strip():
        seqs[nm].append(ln.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
L = len(next(iter(seqs.values())))
ncod = L // 3

# how many codon columns survive COMPLETE deletion across all 22 (what PAML did)?
keep = 0
for c in range(ncod):
    col = [seqs[k][3 * c:3 * c + 3] for k in seqs]
    if all(all(ch in "TCAG" for ch in x) for x in col):
        keep += 1
print("=== C8  site sets ===")
print(f"alignment codons                                     : {ncod}")
print(f"codons surviving complete deletion over all 22 taxa   : {keep}")
print(f"yn00 log reports                                      : 424 sites "
      f"('420 patterns at 424 / 424 sites')")
print(f"-> yn00 (BOTH its NG86 section A and its YN00 section B) analysed "
      f"{keep} codons.")

ms = list(csv.DictReader((GD / "results" / "pairwise_dnds_matrix.tsv").open(), delimiter="\t"))
nc = [int(r["n_codons_compared"]) for r in ms]
print(f"manuscript pairwise-deletion codons per pair           : "
      f"min {min(nc)}, median {int(median(nc))}, max {max(nc)}")
print(f"-> the manuscript estimator used {min(nc)}-{max(nc)} codons per pair, "
      f"yn00 used a fixed {keep}.")
print("   The three estimators share the same 22 SPECIES but NOT the same SITES.")
print("   scripts/04_parse_compare.py docstring says 'Three estimators on the SAME")
print("   22 species' (true) but results/R3.2_estimator_comparison.md says 'Same 22")
print("   species, same codon alignment' (not true at site level).")

# consequence: the two 'undefined' counts in the same .md are about different data
jcundef = sum(1 for r in ms if float(r["ps"]) >= 0.75)
print(f"\n   consequence -- two different 'undefined' counts appear in the SAME .md:")
print(f"     'NG86 undefined 149/231'  <- yn00, {keep}-codon complete-deletion set")
print(f"     'ps >= 0.75: 204/231'     <- manuscript, pairwise-deletion sets "
      f"(recomputed: {jcundef}/231)")
print(f"   Both are individually correct but they are NOT the same quantity on the")
print(f"   same data, and the .md presents them adjacently without saying so.")

# ======================================================= C9 dropout bias
print("\n=== C9  is the NG86 'ordering preserved' result biased by dropout? ===")
tsv = list(csv.DictReader((A / "results" / "pairwise_three_estimators.tsv").open(), delimiter="\t"))
def f(x): return None if x in ("", "None") else float(x)
GROUPS = ["pathogenic-pathogenic", "cross", "environmental-environmental"]
print(f"{'group':32s} {'kept':>7s} {'ps|kept':>9s} {'ps|dropped':>11s} {'yn_dS|kept':>11s} {'yn_dS|drop':>11s}")
for g in GROUPS:
    sub = [r for r in tsv if r["comparison_type"] == g]
    kept = [r for r in sub if f(r["ng_omega"]) is not None]
    drop = [r for r in sub if f(r["ng_omega"]) is None]
    print(f"{g:32s} {len(kept):3d}/{len(sub):3d} "
          f"{mean([f(r['ms_ps']) for r in kept]):9.4f} "
          f"{mean([f(r['ms_ps']) for r in drop]):11.4f} "
          f"{mean([f(r['yn_dS']) for r in kept]):11.4f} "
          f"{mean([f(r['yn_dS']) for r in drop]):11.4f}")
print("""
Retention is 62% / 24% / 38% -- strongly non-random AND non-monotonic across the
three groups being compared, and the kept pairs are the least-diverged ones in
every group. So the NG86 row of the comparison table does not compare the three
groups; it compares three differently-selected low-divergence subsets. The NG86
'ordering preserved' line is therefore the WEAKEST of the three, not a second
independent confirmation. Only the YN00 row (231/231, no dropout) supports the
ordering claim on complete data.""")

# does the ordering survive if we restrict ALL estimators to the NG86-defined subset?
print("=== C9b  common-subset check: all three estimators on the 82 NG86-defined pairs ===")
sub82 = [r for r in tsv if f(r["ng_omega"]) is not None]
for key, nmk in (("ms_omega", "manuscript pn/ps"), ("ng_omega", "NG86"), ("yn_omega", "YN00")):
    data = [[f(r[key]) for r in sub82 if r["comparison_type"] == g] for g in GROUPS]
    H, p = sps.kruskal(*data)
    mus = [mean(d) for d in data]
    oknm = " < ".join(GROUPS[i].split("-")[0][:4] for i in sorted(range(3), key=lambda i: mus[i]))
    print(f"  {nmk:20s} means {mus[0]:.4f} / {mus[1]:.4f} / {mus[2]:.4f}  "
          f"H={H:.3f} p={p:.3e}  ordering path<cross<env? "
          f"{'YES' if mus[0] < mus[1] < mus[2] else 'NO (' + oknm + ')'}")

# ======================================================= C10 YN00 reliability
print("\n=== C10  are the YN00 numbers themselves in a trustworthy range? ===")
dS = [f(r["yn_dS"]) for r in tsv]
om = [f(r["yn_omega"]) for r in tsv]
print(f"YN00 dS over 231 pairs: min {min(dS):.3f}, median {median(dS):.3f}, max {max(dS):.3f}")
print(f"pairs with dS > 1 (i.e. >1 synonymous substitution/site): "
      f"{sum(1 for x in dS if x > 1)}/231")
print(f"pairs with dS > 2: {sum(1 for x in dS if x > 2)}/231")
print(f"YN00 omega: min {min(om):.4f}, max {max(om):.4f}")
print(f"YN00 omega flagged by PAML as uncomputable (-1/99/999): "
      f"{sum(1 for x in om if x in (-1.0, 99.0, 999.0))}/231")
# SE of dS, where available
sb = Path(A / "results" / "yn00" / "yn00.out").read_text().split(
    "(B) Yang & Nielsen (2000) method", 1)[1]
import re
ses = [(float(m.group(1)), float(m.group(2))) for m in
       re.finditer(r"([\d.]+)\s*\+-\s*([\d.]+)\s*$", sb, re.M)]
big = [(d, s) for d, s in ses if s > d]
print(f"pairs whose YN00 dS has SE LARGER than dS itself: {len(big)}/{len(ses)}")
print("""
dS of 3-4 substitutions/site with SEs that exceed the estimate for most pairs
means YN00 is also extrapolating far past the data. YN00 does not rescue the
pairwise approach; it only removes the hard 'undefined' wall. The honest reading
is that pairwise synonymous divergence is saturated at this depth for ALL three
estimators, which is precisely R3.2's and R4.3's objection, and it is why the
tree-based codeml estimate -- not any pairwise estimator -- has to carry the
conclusion.""")
