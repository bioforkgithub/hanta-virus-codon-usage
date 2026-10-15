#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""C16 node support: which nodes are too weak to rely on?  C17 is the
pathogenic-polyphyly claim robust to collapsing weakly supported branches?"""
import itertools, re, sys
from pathlib import Path

A = Path("/PATH/TO/"
         "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12")
PATH10 = {"HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV", "CHOV", "LANV"}
ROD12 = PATH10 | {"PHV", "TULV"}


def parse(nwk):
    s = re.sub(r"\s+", "", nwk.strip().rstrip(";"))
    p = [0]
    def node():
        if s[p[0]] == "(":
            p[0] += 1
            kids = [node()]
            while s[p[0]] == ",":
                p[0] += 1; kids.append(node())
            p[0] += 1
            lbl = ""
            while p[0] < len(s) and s[p[0]] not in "(),":
                lbl += s[p[0]]; p[0] += 1
            return {"kids": kids, "lbl": lbl}
        lbl = ""
        while p[0] < len(s) and s[p[0]] not in "(),":
            lbl += s[p[0]]; p[0] += 1
        return {"kids": [], "lbl": lbl}
    return node()


def leaves(n):
    if not n["kids"]:
        return {n["lbl"].split(":")[0]}
    out = set()
    for k in n["kids"]:
        out |= leaves(k)
    return out


T = parse((A / "results" / "S22_ml.treefile").read_text())
ALL = leaves(T)
assert len(ALL) == 22

print("=== C16  internal-branch support (SH-aLRT / UFBoot) ===")
rows = []
def visit(n):
    if n["kids"]:
        m = re.match(r"([\d.]+)/([\d.]+)", n["lbl"])
        lv = leaves(n)
        if m and 1 < len(lv) < 22:
            rows.append((float(m.group(1)), float(m.group(2)), sorted(lv)))
    for k in n["kids"]:
        visit(k)
visit(T)
rows.sort()
WEAK_A, WEAK_B = 80.0, 95.0   # the usual IQ-TREE thresholds
print(f"{'SH-aLRT':>8s} {'UFBoot':>7s}  verdict        clade")
for a, b, lv in rows:
    strong = a >= WEAK_A and b >= WEAK_B
    v = "strong" if strong else ("WEAK" if (a < WEAK_A and b < WEAK_B) else "marginal")
    print(f"{a:8.1f} {b:7.1f}  {v:13s}  {{{','.join(lv)}}}"
          + ("" if len(lv) <= 8 else f" [{len(lv)} taxa]"))
nweak = sum(1 for a, b, _ in rows if not (a >= WEAK_A and b >= WEAK_B))
print(f"\n{len(rows)} internal branches; {nweak} fail the conventional "
      f"SH-aLRT>=80 AND UFBoot>=95 bar.")

print("""
Branches that matter for the claims actually being made:""")
for a, b, lv in rows:
    s = set(lv)
    tag = None
    if s == {"PHV", "PUUV", "TULV"}:
        tag = ("carries the pathogenic-POLYPHYLY claim: it places the two "
               "rodent-borne environmental species PHV+TULV together with the "
               "pathogenic PUUV")
    elif "RKPV" in s and len(s) == 7:
        tag = ("carries the 'RKPV nests inside the rodent clade' claim")
    elif s == {"PHV", "TULV"}:
        tag = "PHV+TULV sister pair"
    elif s == {"ANDV", "LANV"}:
        tag = "inside the New World HPS clade"
    if tag:
        print(f"  {a:6.1f}/{b:5.1f}  {tag}")

# ---------- C17 robustness of the polyphyly claim
print("\n=== C17  can the 10 pathogenic species be monophyletic under ANY "
      "resolution of the weak branches? ===")

def collapse(n, thr_a, thr_b):
    """return a tree with weakly supported internal edges collapsed into the parent"""
    if not n["kids"]:
        return {"kids": [], "lbl": n["lbl"]}
    kids = []
    for k in n["kids"]:
        ck = collapse(k, thr_a, thr_b)
        if ck["kids"]:
            m = re.match(r"([\d.]+)/([\d.]+)", k["lbl"])
            weak = m is not None and (float(m.group(1)) < thr_a or float(m.group(2)) < thr_b)
            if weak:
                kids.extend(ck["kids"])        # dissolve this node
                continue
        kids.append(ck)
    return {"kids": kids, "lbl": n["lbl"]}


def clade_possible(T, X):
    """Could X be a clade in SOME binary resolution of (possibly multifurcating) T?
    X is realisable iff at some node v, X is exactly the union of the leaf sets of a
    subset of v's children (X may also equal leaves(v) itself)."""
    ok = [False]
    allv = leaves(T)
    def rc(v):
        lv = leaves(v)
        if X <= lv:
            if X == lv:
                ok[0] = True
            elif v["kids"]:
                ksets = [leaves(k) for k in v["kids"]]
                # X must be a union of whole child leaf-sets
                if all((s <= X) or (s & X == set()) for s in ksets):
                    ok[0] = True
            for k in v["kids"]:
                rc(k)
    rc(T)
    # also the complement, since an unrooted bipartition splits either way
    Y = allv - X
    def rc2(v):
        lv = leaves(v)
        if Y <= lv:
            if Y == lv:
                ok[0] = True
            elif v["kids"]:
                ksets = [leaves(k) for k in v["kids"]]
                if all((s <= Y) or (s & Y == set()) for s in ksets):
                    ok[0] = True
            for k in v["kids"]:
                rc2(k)
    rc2(T)
    return ok[0]


for thr in [(0, 0), (80, 95), (0, 95), (0, 70)]:
    C = collapse(T, thr[0], thr[1]) if thr != (0, 0) else T
    ninternal = 0
    def cnt(v):
        global ninternal
        if v["kids"] and 1 < len(leaves(v)) < 22:
            ninternal += 1
        for k in v["kids"]:
            cnt(k)
    ninternal = 0; cnt(C)
    pp = clade_possible(C, PATH10)
    rr = clade_possible(C, ROD12)
    lbl = "ML tree as-is" if thr == (0, 0) else f"collapse SH-aLRT<{thr[0]} or UFBoot<{thr[1]}"
    print(f"  {lbl:36s} internal splits left {ninternal:2d}  "
          f"pathogenic-10 clade possible? {pp}   rodent-12 clade possible? {rr}")

print("""
Reading: even after dissolving every branch that fails UFBoot>=95, the 10
pathogenic species still cannot be made monophyletic, and neither can the 12
rodent-borne ones. The polyphyly facts this analysis rests on are therefore NOT an
artefact of a few weak nodes -- they survive collapsing them. (The specific
'RKPV nests inside the rodent clade' detail does not: see its support above.)""")

# ---------- consensus vs ML tree
con = parse((A / "results" / "S22_ml.contree").read_text())
def bips(T):
    allt = leaves(T); out = set()
    def rc(v):
        s = leaves(v)
        if 1 < len(s) < len(allt):
            a, b = frozenset(s), frozenset(allt - s)
            out.add(a if min(a) < min(b) else b)
        for k in v["kids"]:
            rc(k)
    rc(T); return out
bm, bc = bips(T), bips(con)
print(f"\n=== C18  ML tree vs UFBoot consensus tree ===")
print(f"  ML splits {len(bm)}, consensus splits {len(bc)}, "
      f"symmetric difference {len(bm ^ bc)}")
for s in sorted(bm - bc, key=len):
    print(f"    in ML only: {{{','.join(sorted(s))}}}")
for s in sorted(bc - bm, key=len):
    print(f"    in consensus only: {{{','.join(sorted(s))}}}")
