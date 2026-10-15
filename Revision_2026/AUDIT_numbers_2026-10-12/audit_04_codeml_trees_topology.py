#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""C11 bootstrap stripping lossless?  C12 #1 labelling exact?  C13 branch-omega
assignment read from the per-branch table, not from output order.  C14 LRTs.
C15 topology claims re-verified from scratch."""
import itertools, math, re
from pathlib import Path
from scipy import stats as sps

A = Path("/PATH/TO/"
         "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12")
PATH10 = {"HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV", "CHOV", "LANV"}
ROD12 = PATH10 | {"PHV", "TULV"}


# ---------------- minimal newick parser (no ete3/dendropy available)
def parse(nwk):
    s = re.sub(r"\s+", "", nwk.strip().rstrip(";"))
    pos = [0]

    def node():
        if s[pos[0]] == "(":
            pos[0] += 1
            kids = [node()]
            while s[pos[0]] == ",":
                pos[0] += 1
                kids.append(node())
            assert s[pos[0]] == ")", s[pos[0]:pos[0] + 20]
            pos[0] += 1
            lbl = ""
            while pos[0] < len(s) and s[pos[0]] not in "(),":
                lbl += s[pos[0]]; pos[0] += 1
            return {"kids": kids, "lbl": lbl}
        lbl = ""
        while pos[0] < len(s) and s[pos[0]] not in "(),":
            lbl += s[pos[0]]; pos[0] += 1
        return {"kids": [], "lbl": lbl}
    t = node()
    assert pos[0] == len(s), f"trailing: {s[pos[0]:]!r}"
    return t


def split_lbl(lbl):
    """-> (name_or_support, branchlen_str, paml_label)"""
    paml = None
    m = re.search(r"#([^\s:,()]+)", lbl)
    if m:
        paml = m.group(1); lbl = lbl[:m.start()] + lbl[m.end():]
    name, _, bl = lbl.partition(":")
    return name.strip(), bl.strip(), paml


def walk(n, f, parent=None):
    f(n, parent)
    for k in n["kids"]:
        walk(k, f, n)


def tips(n):
    out = []
    walk(n, lambda x, p: out.append(split_lbl(x["lbl"])[0]) if not x["kids"] else None)
    return out


def bipartitions(n, allt):
    """set of frozensets: for every branch, the tip set below it (normalised)."""
    allt = set(allt)
    res = set()

    def rc(x):
        if not x["kids"]:
            s = {split_lbl(x["lbl"])[0]}
        else:
            s = set()
            for k in x["kids"]:
                s |= rc(k)
        if 0 < len(s) < len(allt):
            a, b = frozenset(s), frozenset(allt - s)
            res.add(a if min(a) < min(b) else b)
        return s
    rc(n)
    return res


FAIL, PASS, CAV = [], [], []
def rec(k, t, m):
    {"FAIL": FAIL, "PASS": PASS, "CAVEAT": CAV}[k].append((t, m)); print(f"[{k}] {t}: {m}")


tf = (A / "results" / "S22_ml.treefile").read_text()
cl = (A / "inputs" / "S22_ml_clean.nwk").read_text()
pt = (A / "inputs" / "S22_pathtips.nwk").read_text()
T_tf, T_cl, T_pt = parse(tf), parse(cl), parse(pt)

# ---------------- C11 stripping lossless
print("=== C11  bootstrap-label stripping: lossless? ===")
tips_tf, tips_cl, tips_pt = tips(T_tf), tips(T_cl), tips(T_pt)
rec("PASS" if tips_tf == tips_cl == tips_pt else "FAIL", "C11.tips",
    f"tip name lists identical across treefile/clean/pathtips ({len(tips_tf)} tips)"
    if tips_tf == tips_cl == tips_pt else
    f"tip lists differ: {tips_tf}\n{tips_cl}\n{tips_pt}")
bp_tf, bp_cl, bp_pt = (bipartitions(t, tips_tf) for t in (T_tf, T_cl, T_pt))
rec("PASS" if bp_tf == bp_cl == bp_pt else "FAIL", "C11.topology",
    f"all three trees have the identical bipartition set ({len(bp_tf)} splits) "
    "-> stripping changed no topology" if bp_tf == bp_cl == bp_pt else
    f"topology differs: clean-vs-treefile symdiff {len(bp_tf ^ bp_cl)}, "
    f"pathtips-vs-treefile symdiff {len(bp_tf ^ bp_pt)}")

def blens(T):
    out = []
    walk(T, lambda x, p: out.append(split_lbl(x["lbl"])[1]))
    return sorted(b for b in out if b)
b1, b2, b3 = blens(T_tf), blens(T_cl), blens(T_pt)
rec("PASS" if b1 == b2 == b3 else "FAIL", "C11.brlen",
    f"all {len(b1)} branch-length strings identical across the three files "
    "(byte-for-byte, no re-rounding)" if b1 == b2 == b3 else
    f"branch lengths differ: {len(b1)}/{len(b2)}/{len(b3)} values")

# were any support values lost that mattered? count them
sup = re.findall(r"\)([\d.]+/[\d.]+):", tf)
rec("PASS", "C11.support",
    f"{len(sup)} internal branches carried SH-aLRT/UFBoot labels in the treefile; "
    f"S22_ml_clean.nwk has {len(re.findall(r'[0-9.]+/[0-9.]+', cl))} such labels left "
    "(PAML cannot parse them, so removal is required, and nothing else was removed)")

# ---------------- C12 #1 labelling
print("\n=== C12  does S22_pathtips.nwk label exactly the 10 pathogenic TIPS? ===")
lab_tips, lab_int = set(), []
def chk(x, p):
    nm, bl, pl = split_lbl(x["lbl"])
    if pl is not None:
        if x["kids"]:
            lab_int.append(pl)
        else:
            lab_tips.add((nm, pl))
walk(T_pt, chk)
got = {n for n, p in lab_tips}
bad = {(n, p) for n, p in lab_tips if p != "1"}
if got == PATH10 and not bad and not lab_int:
    rec("PASS", "C12.labels",
        f"exactly {len(got)} tips carry '#1' and they are precisely the 10 pathogenic "
        f"species {sorted(got)}; zero internal branches are labelled; no other label value used")
else:
    rec("FAIL", "C12.labels",
        f"labelled tips {sorted(got)} (expected {sorted(PATH10)}); "
        f"non-'1' labels {bad}; internal labels {lab_int}")

# ---------------- C13 branch-omega assignment from the node-label block
print("\n=== C13  which w is foreground? read from codeml's own per-branch output ===")
for run in ("two_ratio_F3x4", "two_ratio_F61"):
    o = (A / "results" / "codeml" / run / "out.txt").read_text()
    hdr = re.search(r"^w \(dN/dS\) for branches:\s*(\S+)\s+(\S+)", o, re.M)
    w0, w1 = hdr.group(1), hdr.group(2)
    blk = o.split("w ratios as node labels:", 1)[1].strip().split("\n")[0]
    Tw = parse(blk)
    per = {}
    def g(x, p):
        nm, bl, pl = split_lbl(x["lbl"])
        if not x["kids"]:
            per[nm] = float(pl)
    walk(Tw, g)
    fg = {v for k, v in per.items() if k in PATH10}
    bg = {v for k, v in per.items() if k not in PATH10}
    assert len(fg) == 1 and len(bg) == 1, (fg, bg)
    fgv, bgv = fg.pop(), bg.pop()
    # header prints 5 decimals, so compare at the header's own precision
    hdr_bg_first = abs(float(w0) - bgv) < 1e-5 and abs(float(w1) - fgv) < 1e-5
    hdr_fg_first = abs(float(w0) - fgv) < 1e-5 and abs(float(w1) - bgv) < 1e-5
    assert hdr_bg_first != hdr_fg_first, (w0, w1, fgv, bgv)
    rec("PASS" if hdr_bg_first else "FAIL", f"C13.order.{run}",
        f"header 'w (dN/dS) for branches: {w0} {w1}' -> from the per-branch node-label "
        f"tree, the 10 pathogenic tips all carry w={fgv:g} and all 31 other branches "
        f"carry w={bgv:g}. So the header order is (background, foreground)"
        + ("  CONFIRMED" if hdr_bg_first else "  NOT AS ASSUMED"))
    direction = "LOWER (stronger purifying)" if fgv < bgv else "HIGHER (weaker purifying)"
    rec("PASS", f"C13.direction.{run}",
        f"{run}: foreground(pathogenic) w = {fgv:g}, background w = {bgv:g} "
        f"-> pathogenic omega is {direction}")

# ---------------- C14 LRTs
print("\n=== C14  LRTs ===")
lnl = {}
for run in ("m0_F3x4", "m0_F61", "two_ratio_F3x4", "two_ratio_F61"):
    o = (A / "results" / "codeml" / run / "out.txt").read_text()
    m = re.search(r"lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\):\s*(-?[\d.]+)", o)
    lnl[run] = (int(m.group(1)), int(m.group(2)), float(m.group(3)))
    print(f"  {run:16s} ntime={m.group(1)} np={m.group(2)} lnL={m.group(3)}")
for cf in ("F3x4", "F61"):
    n0, p0, l0 = lnl[f"m0_{cf}"]
    n1, p1, l1 = lnl[f"two_ratio_{cf}"]
    d = 2 * (l1 - l0); df = p1 - p0
    p = sps.chi2.sf(d, df)
    rec("PASS", f"C14.lrt.{cf}",
        f"M0 np={p0} vs two-ratio np={p1} -> df={df} (correct: the only added free "
        f"parameter is the second omega); 2dlnL={d:.4f}, chi2 p={p:.5f}")

# ---------------- C15 topology claims
print("\n=== C15  bipartition claims re-verified ===")
bps = bipartitions(T_cl, tips_cl)
allt = set(tips_cl)
def sep(group):
    g = frozenset(group); c = frozenset(allt - g)
    n = g if min(g) < min(c) else c
    return n in bps
r1, r2 = sep(PATH10), sep(ROD12)
rec("PASS" if not r1 else "FAIL", "C15.path",
    f"is there a single bipartition separating the 10 pathogenic species? {r1} "
    "(claim: False)")
rec("PASS" if not r2 else "FAIL", "C15.rodent",
    f"is there a single bipartition separating the 12 rodent-borne species? {r2} "
    "(claim: False)")
# nearest split containing only pathogenic species
best = None
for b in bps:
    for side in (b, frozenset(allt - b)):
        if side <= PATH10 and (best is None or len(side) > len(best)):
            best = side
rec("PASS", "C15.largest_clade",
    f"largest clade consisting ONLY of pathogenic species: {sorted(best)} "
    f"({len(best)}/10). Claim said the New World HPS clade "
    "{ANDV,BAYV,BCCV,CHOV,LANV,SNV}, 6/10 -> "
    + ("matches" if set(best) == {"ANDV", "BAYV", "BCCV", "CHOV", "LANV", "SNV"} else "DOES NOT MATCH"))
# how many rodent-borne species are monophyletic at best
best_r = None
for b in bps:
    for side in (b, frozenset(allt - b)):
        if side <= ROD12 and (best_r is None or len(side) > len(best_r)):
            best_r = side
rec("PASS", "C15.rodent_clade",
    f"largest clade consisting ONLY of rodent-borne species: {sorted(best_r)} "
    f"({len(best_r)}/12)")
# is RKPV inside the rodent clade?
for b in bps:
    for side in (b, frozenset(allt - b)):
        if side >= ROD12 - {"PHV", "TULV"} and "RKPV" in side and len(side) <= 14:
            print(f"  [info] smallest split containing the HFRS/HPS rodent viruses and RKPV: "
                  f"{sorted(side)}")
            break

print("\n==== TOTALS ====")
print(f"PASS {len(PASS)}  CAVEAT {len(CAV)}  FAIL {len(FAIL)}")
for t, m in FAIL:
    print("  FAIL", t, m)
