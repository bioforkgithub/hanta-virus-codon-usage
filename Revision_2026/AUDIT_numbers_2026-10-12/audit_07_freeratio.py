#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""C20: free-ratio F3x4 finished. Does the per-branch omega distribution support the
two-ratio result, and is the two-ratio model itself an adequate description?"""
import re
from pathlib import Path
from statistics import mean, median
from scipy import stats as sps

A = Path("/PATH/TO/"
         "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12")
PATH10 = {"HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV", "CHOV", "LANV"}

o = (A / "results" / "codeml" / "free_ratio_F3x4" / "out.txt").read_text()
## the node-label tree can WRAP over several lines: take everything up to the ';'
tail = o.split("w ratios as node labels:", 1)[1]
blk = tail[:tail.index(";") + 1]

# pull every "NAME #value" for tip names. Match against the ACTUAL taxon list:
# PHV and SNV are 3 letters while every other abbreviation is 4, so a [A-Z]{4}
# pattern silently drops exactly those two.
import csv
TAXA = [r["abbrev"] for r in csv.DictReader(
    (A / "inputs" / "labels.tsv").open(), delimiter="\t")]
assert len(TAXA) == 22
per = dict((m.group(1), float(m.group(2))) for m in re.finditer(
    r"\b(" + "|".join(sorted(TAXA, key=len, reverse=True)) + r")\s*#([\d.eE+-]+)", blk))
print(f"terminal branches recovered: {len(per)} of 22")
missing = set(TAXA) - set(per)
assert not missing, f"FAILED TO PARSE: {missing}"
fg = {k: v for k, v in per.items() if k in PATH10}
bg = {k: v for k, v in per.items() if k not in PATH10}
print(f"\npathogenic terminals (n={len(fg)}):")
for k, v in sorted(fg.items(), key=lambda kv: kv[1]):
    print(f"   {k}  {v:.5f}")
print(f"\nenvironmental terminals (n={len(bg)}):")
for k, v in sorted(bg.items(), key=lambda kv: kv[1]):
    print(f"   {k}  {v:.5f}")
a, b = list(fg.values()), list(bg.values())
print(f"\npathogenic   terminal w: mean {mean(a):.5f}  median {median(a):.5f}")
print(f"environmental terminal w: mean {mean(b):.5f}  median {median(b):.5f}")
U, pu = sps.mannwhitneyu(a, b, alternative="two-sided")
print(f"Mann-Whitney U={U:.1f}  p={pu:.4f}  (two-sided, 10 vs 12 terminal branches)")
print(f"direction: pathogenic terminals have "
      f"{'LOWER' if mean(a) < mean(b) else 'HIGHER'} mean omega")

# --------------------------------------------------- is the free-ratio fit usable?
sec = o.split("dN & dS for each branch", 1)[1].split("tree length for dN", 1)[0]
# columns: branch  t  N  S  dN/dS  dN  dS  N*dN  S*dS
br = [dict(a=int(m.group(1)), b=int(m.group(2)), t=float(m.group(3)),
           w=float(m.group(6)), dS=float(m.group(8)))
      for m in re.finditer(
          r"^\s*(\d+)\.\.(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)",
          sec, re.M)]
assert len(br) == 41, len(br)
term = [r for r in br if r["b"] <= 22]
intr = [r for r in br if r["b"] > 22]
print("\n--- is the free-ratio fit usable? ---")
for nm, g in (("terminal", term), ("internal", intr)):
    print(f"  {nm:8s} n={len(g):2d}  mean w={mean(r['w'] for r in g):10.5f}  "
          f"median w={median(r['w'] for r in g):9.5f}  "
          f"mean dS={mean(r['dS'] for r in g):7.3f}")
deg = [r for r in br if r["dS"] < 0.01]
print(f"  branches with dS < 0.01, where omega is UNIDENTIFIABLE: {len(deg)}/41; "
      f"all internal = {all(r['b'] > 22 for r in deg)}")
print(f"  their omegas: " + ", ".join(f"{r['w']:.1f}" for r in deg))
print("""
  => the free-ratio fit is DEGENERATE on internal branches: with 41 free omegas on
     448 saturated codons the optimiser drives internal dS to ~0 and omega explodes
     (up to 289.7). Internal free-ratio omegas must NOT be reported, and any model
     comparison AGAINST the free-ratio model is uninformative for the same reason.
     The 22 TERMINAL branches are well behaved (dS 2-30, omega 0.001-0.067) and the
     pathogenic-vs-environmental comparison above uses only those.""")

# nested LRTs
lnl = {}
for run in ("m0_F3x4", "two_ratio_F3x4", "free_ratio_F3x4"):
    t = (A / "results" / "codeml" / run / "out.txt").read_text()
    m = re.search(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?[\d.]+)", t)
    lnl[run] = (int(m.group(1)), float(m.group(2)))
print("\n--- nested LRTs (all CodonFreq=2 / F3x4, same tree) ---")
for n0, n1 in (("m0_F3x4", "two_ratio_F3x4"), ("m0_F3x4", "free_ratio_F3x4"),
               ("two_ratio_F3x4", "free_ratio_F3x4")):
    p0, l0 = lnl[n0]; p1, l1 = lnl[n1]
    d = 2 * (l1 - l0); df = p1 - p0
    note = "  <- UNINFORMATIVE (free-ratio fit is degenerate)" if "free" in n1 else ""
    print(f"  {n0:16s} -> {n1:16s} 2dlnL={d:9.3f}  df={df:3d}  "
          f"p={sps.chi2.sf(d, df):.3e}{note}")
