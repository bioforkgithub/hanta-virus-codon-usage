#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""C21: independently verify the orchestrator's four-model CodonFreq sweep and its
degeneracy diagnostics, in results/CRITICAL_codonfreq_direction_flip.md."""
import csv, re
from pathlib import Path
from statistics import mean, median
from scipy import stats as sps

A = Path("/PATH/TO/"
         "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12")
PATH10 = {"HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV", "CHOV", "LANV"}
TAXA = [r["abbrev"] for r in csv.DictReader(
    (A / "inputs" / "labels.tsv").open(), delimiter="\t")]
NAMEPAT = "|".join(sorted(TAXA, key=len, reverse=True))
TBOUND = 50.0     # codeml's documented upper bound on branch length t


def branches(out):
    sec = out.split("dN & dS for each branch", 1)[1].split("tree length for dN", 1)[0]
    return [dict(a=int(m.group(1)), b=int(m.group(2)), t=float(m.group(3)),
                 w=float(m.group(6)), dN=float(m.group(7)), dS=float(m.group(8)))
            for m in re.finditer(
                r"^\s*(\d+)\.\.(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)",
                sec, re.M)]


def fgbg(out):
    """foreground/background w read from the per-branch node-label tree, by taxon name"""
    tail = out.split("w ratios as node labels:", 1)[1]
    blk = tail[:tail.index(";") + 1]
    per = dict((m.group(1), float(m.group(2))) for m in re.finditer(
        r"\b(" + NAMEPAT + r")\s*#([\d.eE+-]+)", blk))
    assert len(per) == 22, sorted(set(TAXA) - set(per))
    fg = {v for k, v in per.items() if k in PATH10}
    bg = {v for k, v in per.items() if k not in PATH10}
    return (fg.pop() if len(fg) == 1 else None), (bg.pop() if len(bg) == 1 else None), per


def lnl_np(out):
    m = re.search(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?[\d.]+)", out)
    return int(m.group(1)), float(m.group(2))


print("=== C21a  the four-model sweep, re-read from raw out.txt ===")
print(f"{'CodonFreq':>10s} {'M0 lnL':>13s} {'2ratio lnL':>13s} {'w bg':>9s} {'w fg':>9s} "
      f"{'dir':>9s} {'treelen':>9s} {'dS len':>9s} {'t@bound':>8s} {'maxdS':>8s}")
rows = {}
sweep = A / "results" / "codeml_sweep"
for cf in range(4):
    cands = [sweep / f"two_ratio_cf{cf}", sweep / f"m0_cf{cf}"]
    tr = cands[0] / "out.txt"
    m0 = cands[1] / "out.txt"
    if not tr.exists():
        # fall back to the original four runs for cf 2 and 3
        alt = {2: ("two_ratio_F3x4", "m0_F3x4"), 3: ("two_ratio_F61", "m0_F61")}
        if cf in alt:
            tr = A / "results" / "codeml" / alt[cf][0] / "out.txt"
            m0 = A / "results" / "codeml" / alt[cf][1] / "out.txt"
    if not tr.exists():
        print(f"{cf:>10d}   MISSING two-ratio run")
        continue
    o = tr.read_text()
    fg, bg, per = fgbg(o)
    br = branches(o)
    tl = sum(r["t"] for r in br)
    dsl = sum(r["dS"] for r in br)
    nb = sum(1 for r in br if r["t"] >= TBOUND - 1e-6)
    _, l1 = lnl_np(o)
    l0 = lnl_np(m0.read_text())[1] if m0.exists() else float("nan")
    rows[cf] = dict(l0=l0, l1=l1, fg=fg, bg=bg, tl=tl, dsl=dsl, nb=nb,
                    maxdS=max(r["dS"] for r in br), per=per)
    print(f"{cf:>10d} {l0:13.2f} {l1:13.2f} {bg:9.5f} {fg:9.5f} "
          f"{'LOWER' if fg < bg else 'HIGHER':>9s} {tl:9.1f} {dsl:9.1f} {nb:8d} "
          f"{max(r['dS'] for r in br):8.2f}")

print("""
'LOWER' = pathogenic omega below background = supports the manuscript's direction.
't@bound' counts branches pinned at codeml's t upper bound of 50.""")

print("\n=== C21b  LRTs for each CodonFreq ===")
for cf, d in sorted(rows.items()):
    dl = 2 * (d["l1"] - d["l0"])
    print(f"  CodonFreq={cf}  2dlnL={dl:7.3f}  df=1  p={sps.chi2.sf(dl, 1):.5f}")

print("\n=== C21c  free_ratio_F3x4 degeneracy diagnostics (their numbers vs mine) ===")
o = (A / "results" / "codeml" / "free_ratio_F3x4" / "out.txt").read_text()
br = branches(o)
tl = sum(r["t"] for r in br)
nb = sum(1 for r in br if r["t"] >= TBOUND - 1e-6)
ws = [r["w"] for r in br]
print(f"  branches            : {len(br)} (expect 41)")
print(f"  tree length (sum t) : {tl:.1f}      [they report 260.4]")
print(f"  branches at t bound : {nb}          [they report 4]")
print(f"  max dS              : {max(r['dS'] for r in br):.1f}       [they report 69.0]")
print(f"  omega range         : {min(ws):.5f} to {max(ws):.1f}   [they report 0.0125 to 289.7]")
lo = sorted(br, key=lambda r: r["w"])[:4]
print("  the 4 LOWEST omegas and whether terminal:")
for r in lo:
    nm = TAXA[r["b"] - 1] if r["b"] <= 22 else f"node{r['b']}"
    print(f"     {r['a']}..{r['b']:<3d} {nm:6s} w={r['w']:.5f}  t={r['t']:.3f}  dS={r['dS']:.3f}"
          f"  {'TERMINAL' if r['b'] <= 22 else 'internal'}")
print("""
  => the omega MINIMUM is the point of disagreement. Checked below against the
     node-label tree, which labels by taxon name and cannot be mis-indexed.""")
fg, bg, per = fgbg(o)
print(f"  node-label tree, 4 lowest terminal omegas: "
      + ", ".join(f"{k} {v:.5f}" for k, v in sorted(per.items(), key=lambda kv: kv[1])[:4]))
print(f"  per-branch table agrees with node-label tree for all 22 terminals: "
      f"{all(abs(per[TAXA[r['b']-1]] - r['w']) < 2e-4 for r in br if r['b'] <= 22)}")
