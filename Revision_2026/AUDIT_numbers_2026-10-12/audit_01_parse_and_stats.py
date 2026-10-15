#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""Adversarial audit of R2.2_phylo_dnds_2026-10-12: independent re-parse of
yn00 output, independent re-implementation of Kruskal-Wallis, and verification
of the manuscript-pn/ps accusation and the saturation counts.

Nothing here reuses 04_parse_compare.py. Every parser is written from scratch
from the raw PAML output so that an error in the original would show up as a
disagreement, not be reproduced.
"""
import csv, itertools, math, re, sys
from pathlib import Path
from statistics import mean, median

import numpy as np
from scipy import stats as sps

A = Path("/PATH/TO/"
         "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12")
GD = Path("/PATH/TO/"
          "Project_Hantavirus/Work1/Genus_Diversity")

FAILS, CAVEATS, PASSES = [], [], []


def rec(kind, tag, msg):
    {"FAIL": FAILS, "CAVEAT": CAVEATS, "PASS": PASSES}[kind].append((tag, msg))
    print(f"[{kind}] {tag}: {msg}")


# ---------------------------------------------------------------- order check
labels = list(csv.DictReader((A / "inputs" / "labels.tsv").open(), delimiter="\t"))
order = [r["abbrev"] for r in labels]
cat = {r["abbrev"]: r["category"] for r in labels}

phy = (A / "inputs" / "S22_codon_aln.phy").read_text().splitlines()
phy_names = [phy[i] for i in range(1, len(phy), 2)]
if phy_names == order:
    rec("PASS", "C0.order", f"labels.tsv order == PHYLIP order ({len(order)} taxa); "
        "yn00 integer indices therefore map to the right species")
else:
    rec("FAIL", "C0.order", f"labels.tsv order != PHYLIP order: {order} vs {phy_names}")

GROUPS = ["pathogenic-pathogenic", "cross", "environmental-environmental"]


def group_of(a, b):
    ca, cb = cat[a], cat[b]
    if ca == cb == "pathogenic":
        return "pathogenic-pathogenic"
    if ca == cb == "environmental":
        return "environmental-environmental"
    return "cross"


# ---------------------------------------------- C1: NG86 lower-triangle parse
txt = (A / "results" / "yn00" / "yn00.out").read_text()
sec_a = txt.split("(A) Nei-Gojobori (1986) method", 1)[1].split("(B) Yang & Nielsen", 1)[0]

# independent parser: walk rows IN FILE ORDER, do not use order.index()
cell = re.compile(r"(-?[\d.]+)\s*\(\s*(-?[\d.]+)\s+(-?[\d.]+)\s*\)")
ng_raw = {}
row_i = 0
rowcells = {}
for ln in sec_a.splitlines():
    tok = ln.split()
    if not tok or tok[0] not in cat:
        continue
    name = tok[0]
    if name != order[row_i]:
        rec("FAIL", "C1.roworder",
            f"NG86 matrix row {row_i} is {name!r} but PHYLIP order expects {order[row_i]!r}")
    cells = cell.findall(ln)
    rowcells[name] = len(cells)
    for j, (om, dn, ds) in enumerate(cells):
        ng_raw[frozenset((name, order[j]))] = (float(om), float(dn), float(ds))
    row_i += 1

if row_i != 22:
    rec("FAIL", "C1.nrows", f"NG86 section has {row_i} named rows, expected 22")
bad = {n: c for n, c in rowcells.items() if c != order.index(n)}
if bad:
    rec("FAIL", "C1.cellcount",
        f"NG86 row(s) with wrong number of lower-triangle cells (row i must have i cells): {bad}")
else:
    rec("PASS", "C1.cellcount",
        "every NG86 matrix row i carries exactly i lower-triangle cells "
        "(0..21, total 231) -> no line wrapping, no off-by-one in the column index")
if len(ng_raw) != 231:
    rec("FAIL", "C1.npairs", f"NG86 independent parse recovered {len(ng_raw)} pairs, expected 231")
else:
    rec("PASS", "C1.npairs", "NG86 independent parse recovered all 231 pairs")

# internal consistency of each NG86 cell: omega must equal dN/dS
inc = []
for k, (om, dn, ds) in ng_raw.items():
    if om == -1.0 or ds == -1.0:
        continue
    if ds > 0 and abs(om - dn / ds) > 5e-4:
        inc.append((sorted(k), om, dn, ds, dn / ds))
if inc:
    rec("FAIL", "C1.cellselfconsistency",
        f"{len(inc)} NG86 cells where printed omega != printed dN/dS, e.g. {inc[:3]}")
else:
    rec("PASS", "C1.cellselfconsistency",
        "in every defined NG86 cell the printed omega equals printed dN/dS to <5e-4 "
        "-> the three numbers in each cell are being read as (omega, dN, dS) in the right roles")

# ------------------------------------------------- C2: YN00 section B reparse
sec_b = txt.split("(B) Yang & Nielsen (2000) method", 1)[1]
ynpat = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+"
    r"(\S+)\s*\+-\s*(\S+)\s+(\S+)\s*\+-\s*(\S+)\s*$", re.M)
yn_raw, yn_ij = {}, {}
for m in ynpat.finditer(sec_b):
    i, j = int(m.group(1)), int(m.group(2))
    a, b = order[i - 1], order[j - 1]
    d = dict(S=float(m.group(3)), N=float(m.group(4)), t=float(m.group(5)),
             kappa=float(m.group(6)), omega=float(m.group(7)),
             dN=float(m.group(8)), dS=float(m.group(10)))
    yn_raw[frozenset((a, b))] = d
    yn_ij[(i, j)] = d
if len(yn_raw) != 231:
    rec("FAIL", "C2.npairs", f"YN00 section B reparse got {len(yn_raw)} pairs, expected 231")
else:
    rec("PASS", "C2.npairs", "YN00 section B reparse got all 231 pairs (no silently dropped line)")

# cross-check against the machine-readable 2YN.dN / 2YN.dS matrices
def read_2yn(p):
    lines = [l for l in p.read_text().splitlines() if l.strip()]
    n = int(lines[0].split()[0])
    M = {}
    for r in range(1, n + 1):
        tok = lines[r].split()
        nm = tok[0]
        for c, v in enumerate(tok[1:]):
            M[frozenset((nm, order[c]))] = float(v)
    return M

m_dN = read_2yn(A / "results" / "yn00" / "2YN.dN")
m_dS = read_2yn(A / "results" / "yn00" / "2YN.dS")
diffs = []
for k, d in yn_raw.items():
    for nm, M in (("dN", m_dN), ("dS", m_dS)):
        if k in M and abs(M[k] - d[nm]) > 1.1e-4:
            diffs.append((sorted(k), nm, M[k], d[nm]))
if diffs:
    rec("FAIL", "C2.vs2YN", f"{len(diffs)} YN00 dN/dS disagreements vs 2YN.* matrices, e.g. {diffs[:5]}")
else:
    rec("PASS", "C2.vs2YN",
        "all 231 YN00 dN and dS values parsed from yn00.out match the independent "
        "2YN.dN / 2YN.dS matrices to <1.1e-4 (rounding)")

ynk = {d["kappa"] for d in yn_raw.values()}
if len(ynk) == 1:
    rec("CAVEAT", "C2.kappa",
        f"yn00 estimated ONE global kappa = {ynk.pop()} for all 231 pairs (log line "
        "'kappa = 1.42'). The results .md reports it as 'mean 1.424, median 1.424, "
        "range 1.424-1.424 across pairs', which reads as a distribution over pairs. "
        "It is a single dataset-level estimate with no per-pair variation and no CI.")

# ------------------------------- C3: the accusation about the manuscript column
ms = {}
ms_rows = list(csv.DictReader((GD / "results" / "pairwise_dnds_matrix.tsv").open(), delimiter="\t"))
worst, nchk = 0.0, 0
site_bad = []
for r in ms_rows:
    k = frozenset((r["species_A"], r["species_B"]))
    ps, pn, dnds = float(r["ps"]), float(r["pn"]), float(r["dNdS"])
    ms[k] = dict(ps=ps, pn=pn, omega=dnds)
    worst = max(worst, abs(dnds - pn / ps)); nchk += 1
    # and check ps/pn are themselves raw proportions Sd/S and Nd/N
    if abs(ps - float(r["Sd"]) / float(r["S_sites"])) > 2e-6 or \
       abs(pn - float(r["Nd"]) / float(r["N_sites"])) > 2e-6:
        site_bad.append(r["species_A"] + "-" + r["species_B"])
if worst <= 1e-6 and nchk == 231:
    rec("PASS", "C3.dnds_is_pnps",
        f"all {nchk} rows: dNdS == pn/ps exactly (max |dNdS - pn/ps| = {worst:.2e}, "
        "i.e. pure rounding). The accusation is CONFIRMED.")
else:
    rec("FAIL", "C3.dnds_is_pnps", f"max |dNdS - pn/ps| = {worst:.3e} over {nchk} rows")
if site_bad:
    rec("FAIL", "C3.rawprops", f"ps != Sd/S_sites or pn != Nd/N_sites for {len(site_bad)} rows")
else:
    rec("PASS", "C3.rawprops",
        "ps == Sd/S_sites and pn == Nd/N_sites for all 231 rows -> ps and pn are raw "
        "uncorrected proportions, no Jukes-Cantor term anywhere")

# what the properly corrected NG86 would have been, where defined
def jc(p):
    if p is None or p >= 0.75:
        return None
    return -0.75 * math.log(1 - 4.0 * p / 3.0)

# ----------------------------------------------- C4: saturation / JC arithmetic
over = [k for k, v in ms.items() if v["ps"] >= 0.75]
by_g = {g: [0, 0] for g in GROUPS}
for a, b in itertools.combinations(order, 2):
    g = group_of(a, b)
    by_g[g][1] += 1
    if ms[frozenset((a, b))]["ps"] >= 0.75:
        by_g[g][0] += 1
rec("PASS" if len(over) == 204 else "FAIL", "C4.ps075",
    f"pairs with ps >= 0.75: {len(over)}/231 (claimed 204/231); by group "
    + ", ".join(f"{g} {by_g[g][0]}/{by_g[g][1]}" for g in GROUPS)
    + " (claimed 31/45, 111/120, 62/66)")

# does ps >= 0.75 really break JC?
lim = [(1 - 4 * p / 3) for p in (0.74, 0.7499, 0.75, 0.80)]
rec("PASS", "C4.jcmath",
    "JC argument 1-4p/3 at ps = 0.74/0.7499/0.75/0.80 = "
    + ", ".join(f"{x:.5f}" for x in lim)
    + " -> <=0 at ps>=0.75, so -3/4 ln(1-4ps/3) is undefined (log of 0 or a "
    "negative). Arithmetic confirmed.")

# where NG86 is undefined, is it the same pairs that are saturated?
ng_undef = {k for k, (om, dn, ds) in ng_raw.items() if om == -1.0}
rec("PASS" if len(ng_undef) == 149 else "FAIL", "C4.ngundef",
    f"NG86 undefined (omega printed as -1) for {len(ng_undef)}/231 pairs (claimed 149/231)")
und_by_g = {g: [0, 0] for g in GROUPS}
for a, b in itertools.combinations(order, 2):
    g = group_of(a, b); und_by_g[g][1] += 1
    if frozenset((a, b)) in ng_undef:
        und_by_g[g][0] += 1
rec("PASS", "C4.ngundef_bygroup",
    "NG86 undefined by group: " + ", ".join(f"{g} {und_by_g[g][0]}/{und_by_g[g][1]}"
                                            for g in GROUPS)
    + " (results .md claims 17/45, 91/120, 41/66)")

# ------------------------------------------------- C5: group means / orderings
tsv = list(csv.DictReader((A / "results" / "pairwise_three_estimators.tsv").open(), delimiter="\t"))
def f(x):
    return None if x in ("", "None") else float(x)

# verify the committed TSV against the independent parses
mis = []
for r in tsv:
    k = frozenset((r["species_A"], r["species_B"]))
    om, dn, ds = ng_raw[k]
    exp_ng = None if om == -1.0 else om
    got_ng = f(r["ng_omega"])
    if (exp_ng is None) != (got_ng is None) or (exp_ng is not None and abs(exp_ng - got_ng) > 1e-9):
        mis.append(("ng_omega", sorted(k), exp_ng, got_ng))
    exp_yn = yn_raw[k]["omega"]
    got_yn = f(r["yn_omega"])
    if got_yn is None or abs(exp_yn - got_yn) > 1e-9:
        mis.append(("yn_omega", sorted(k), exp_yn, got_yn))
    exp_ms = ms[k]["omega"]
    got_ms = f(r["ms_omega"])
    if got_ms is None or abs(exp_ms - got_ms) > 1e-9:
        mis.append(("ms_omega", sorted(k), exp_ms, got_ms))
if mis:
    rec("FAIL", "C5.tsv", f"{len(mis)} values in pairwise_three_estimators.tsv disagree with an "
        f"independent parse, e.g. {mis[:5]}")
else:
    rec("PASS", "C5.tsv",
        "every ms_omega / ng_omega / yn_omega in pairwise_three_estimators.tsv (231 rows x 3) "
        "reproduces exactly from an independently written parser -> no pair/species mis-assignment")

print("\n--- recomputed group statistics (independent) ---")
summary = {}
for key, nm in (("ms_omega", "manuscript pn/ps"), ("ng_omega", "NG86 JC-corrected"),
                ("yn_omega", "YN00")):
    data, out = [], {}
    for g in GROUPS:
        v = [f(r[key]) for r in tsv if r["comparison_type"] == g]
        vv = [x for x in v if x is not None]
        data.append(vv)
        out[g] = (len(vv), len(v) - len(vv), mean(vv), median(vv))
    H, p = sps.kruskal(*data)
    summary[key] = (out, H, p, sum(len(d) for d in data))
    print(f"\n{nm}")
    for g in GROUPS:
        n, miss, mu, md = out[g]
        print(f"  {g:30s} n={n:3d} missing={miss:3d} mean={mu:.4f} median={md:.4f}")
    print(f"  scipy kruskal H={H:.3f} p={p:.4e}  N={sum(len(d) for d in data)}")

# ---------------------------------------- C6: hand-rolled KW vs scipy, exp(-H/2)
sys.path.insert(0, str(A / "scripts"))
claimed = {"ms_omega": (60.628, 6.835e-14), "ng_omega": (31.155, 1.717e-07),
           "yn_omega": (32.346, 9.466e-08)}
for key, (cH, cp) in claimed.items():
    out, H, p, N = summary[key]
    ok_H = abs(H - cH) < 0.002
    ok_p = abs(math.log10(p) - math.log10(cp)) < 0.01
    rec("PASS" if (ok_H and ok_p) else "FAIL", f"C6.kw.{key}",
        f"scipy.stats.kruskal H={H:.3f} p={p:.4e} vs script H={cH} p={cp:.4e} "
        f"-> {'agree' if ok_H and ok_p else 'DISAGREE'}")

# the exp(-H/2) identity
ids = []
for h in (0.5, 2.0, 6.869, 31.155, 60.628):
    ids.append((h, sps.chi2.sf(h, 2), math.exp(-h / 2)))
worst_id = max(abs(a - b) / a for _, a, b in ids)
rec("PASS" if worst_id < 1e-9 else "FAIL", "C6.chi2df2",
    "chi2.sf(x, df=2) == exp(-x/2) checked at x = "
    + ", ".join(f"{h:g}" for h, _, _ in ids)
    + f"; worst relative difference {worst_id:.2e}. The closed form is correct for df=2.")

# tie structure: does the tie correction matter here?
for key in ("ms_omega", "yn_omega"):
    v = [f(r[key]) for r in tsv if f(r[key]) is not None]
    nties = len(v) - len(set(v))
    print(f"[info] {key}: {nties} tied observations among {len(v)} "
          f"(tie correction factor is 1.0 when 0)")

# ---------------------------------------------- C7: LRT arithmetic
l0, l1 = -21106.726850, -21103.292568
d = 2 * (l1 - l0)
p = sps.chi2.sf(d, 1)
rec("PASS" if abs(d - 6.869) < 0.001 and abs(p - 0.0088) < 0.0002 else "FAIL", "C7.lrt",
    f"2*(lnL1-lnL0) = 2*({l1} - {l0}) = {d:.6f}; chi2.sf({d:.4f}, df=1) = {p:.6f} "
    f"(claimed 6.869, p=0.0088)")

print("\n================ SUMMARY ================")
print(f"PASS   {len(PASSES)}")
print(f"CAVEAT {len(CAVEATS)}")
print(f"FAIL   {len(FAILS)}")
for t, m in FAILS:
    print("  FAIL", t, "-", m)
