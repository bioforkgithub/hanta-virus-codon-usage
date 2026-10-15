#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""R2.2 crosstalk fix (2026-10-12, post-audit).

Recompute the free-ratio terminal-omega contrast the way the adversarial audit
(analysis/AUDIT_numbers_2026-10-12.md section E1-bis) requires: EXCLUDE the
terminal branches pinned at codeml's t upper bound (t = 50), where omega is not
an estimate. Then ask how much of the surviving environmental signal is carried
by PHV + TULV, the only two rodent-borne environmental species (R3.1 confound).

Per-branch omegas are parsed out of codeml's own branch table; tip numbers are
mapped to taxon names from the alignment order, and categories come from
inputs/labels.tsv (never from a fixed-width name pattern -- audit hygiene note).
Read-only: nothing is re-run, nothing under Work1/Work2 is touched.
"""
import json, os, re, sys
from scipy.stats import mannwhitneyu

BASE = "/PATH/TO/Project_Hantavirus/Reviews_ans/analysis"
RUN = os.path.join(BASE, "R2.2_phylo_dnds_2026-10-12")
OUT = os.path.join(BASE, "R2.2_interpretation_2026-10-12", "results",
                   "freeratio_boundary_corrected.json")

# ---- taxon order (= codeml tip numbering) from the phylip alignment
order = []
with open(os.path.join(RUN, "inputs", "S22_codon_aln.phy")) as fh:
    fh.readline()
    for line in fh:
        tok = line.split()
        if tok and not re.fullmatch(r"[ACGTNacgtn-]+", tok[0]):
            order.append(tok[0])
assert len(order) == 22, order
tip = {i + 1: name for i, name in enumerate(order)}

# ---- categories
cat = {}
with open(os.path.join(RUN, "inputs", "labels.tsv")) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ia, ic = hdr.index("abbrev"), hdr.index("category")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) > max(ia, ic):
            cat[f[ia]] = f[ic]

# ---- codeml branch table
rows = []
with open(os.path.join(RUN, "results", "codeml", "free_ratio_F3x4", "out.txt")) as fh:
    txt = fh.read()
m = re.search(r"branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS.*?\n(.*?)\n\s*\n", txt, re.S)
for line in m.group(1).strip().splitlines():
    f = line.split()
    if len(f) < 7 or ".." not in f[0]:
        continue
    a, b = (int(x) for x in f[0].split(".."))
    rows.append({"branch": f[0], "child": b, "t": float(f[1]),
                 "omega": float(f[4]), "dS": float(f[6])})
assert len(rows) == 41, len(rows)

term = [r | {"taxon": tip[r["child"]], "category": cat[tip[r["child"]]]}
        for r in rows if r["child"] in tip]
assert len(term) == 22, len(term)

def stat(v):
    s = sorted(v); n = len(s)
    return {"n": n, "mean": round(sum(s) / n, 5),
            "median": round(s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2, 5),
            "min": min(s), "max": max(s)}

boundary = [r for r in term if r["t"] >= 50.0]
path = [r for r in term if r["category"] == "pathogenic"]
env_all = [r for r in term if r["category"] != "pathogenic"]
env_keep = [r for r in env_all if r["t"] < 50.0]
env_shrew = [r for r in env_keep if r["taxon"] not in ("PHV", "TULV")]

res = {"boundary_branches": [{k: r[k] for k in ("taxon", "category", "t", "dS", "omega")}
                             for r in boundary],
       "scenarios": {}}
for name, grp in (("env_all_incl_boundary", env_all),
                  ("env_boundary_dropped", env_keep),
                  ("env_boundary_dropped_minus_PHV_TULV", env_shrew)):
    u, pv = mannwhitneyu([r["omega"] for r in path], [r["omega"] for r in grp],
                         alternative="two-sided")
    res["scenarios"][name] = {
        "pathogenic": stat([r["omega"] for r in path]),
        "environmental": stat([r["omega"] for r in grp]),
        "env_taxa": sorted(r["taxon"] for r in grp),
        "mannwhitney_U": u, "p": round(float(pv), 4)}
res["terminal_dS"] = {
    "pathogenic": stat([r["dS"] for r in path]),
    "environmental_all": stat([r["dS"] for r in env_all]),
    "environmental_boundary_dropped": stat([r["dS"] for r in env_keep])}
res["source"] = "results/codeml/free_ratio_F3x4/out.txt (read-only)"

with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1)
json.dump(res, sys.stdout, indent=1)
print()
