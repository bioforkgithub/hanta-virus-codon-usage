# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""
R4.3 — assemble every codeml run available on 2026-10-12 into one table, so the
branch-model result is reported without selecting a favourite.
Read-only over analysis/R2.2_phylo_dnds_2026-10-12/results/codeml*/.
"""
import re
from pathlib import Path
from scipy.stats import chi2

C = Path("/PATH/TO/"
         "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12/results")
OUT = Path(__file__).resolve().parent.parent / "results"

CF = {"cf0": ("0", "equal, 1/61"), "cf1": ("1", "F1x4"),
      "F3x4": ("2", "F3x4"), "F61": ("3", "F61")}


def parse(p):
    if not p.exists():
        return None
    t = p.read_text()
    d = {}
    m = re.search(r"^lnL\(ntime:\s*(\d+)\s*np:\s*(\d+)\):\s*(-?[\d.]+)", t, re.M)
    if not m:
        return None
    d["ntime"], d["np"], d["lnL"] = int(m.group(1)), int(m.group(2)), float(m.group(3))
    m = re.search(r"^omega \(dN/dS\) =\s*([\d.]+)", t, re.M)
    if m:
        d["omega"] = float(m.group(1))
    m = re.search(r"^w \(dN/dS\) for branches:\s*([\d.]+)\s+([\d.]+)", t, re.M)
    if m:
        d["w_bg"], d["w_fg"] = float(m.group(1)), float(m.group(2))
    m = re.search(r"tree length for dS:\s*([\d.]+)", t)
    if m:
        d["dS_len"] = float(m.group(1))
    m = re.search(r"tree length for dN:\s*([\d.]+)", t)
    if m:
        d["dN_len"] = float(m.group(1))
    # branch table rows look like "  23..1      2.868  1047.7 ..." - take col 2
    ts = []
    seen = False
    for line in t.splitlines():
        if "dN & dS for each branch" in line:
            seen = True
            continue
        if seen:
            f = line.split()
            if len(f) >= 8 and re.fullmatch(r"\d+\.\.\d+", f[0]):
                try:
                    ts.append(float(f[1]))
                except ValueError:
                    pass
    d["max_t"] = max(ts) if ts else None
    return d


runs = {}
for base in (C / "codeml", C / "codeml_sweep"):
    if not base.exists():
        continue
    for sub in sorted(base.iterdir()):
        if sub.is_dir():
            runs[sub.name] = (parse(sub / "out.txt"), sub)

rep = []


def say(s=""):
    print(s)
    rep.append(s)


say("# R4.3 — every codeml branch-model run available, 2026-10-12")
say("")
say("Read-only from `analysis/R2.2_phylo_dnds_2026-10-12/results/codeml/` and "
    "`.../codeml_sweep/`. Nothing re-run, nothing deleted.")
say("")
say("Same ML tree (`S22_ml.treefile`), same alignment (448 codons x 22 species), "
    "same foreground (the 10 pathogenic TERMINAL branches). The ONLY thing that "
    "differs between rows is codeml's `CodonFreq` setting.")
say("")
say("| CodonFreq | M0 lnL | M0 omega | two-ratio lnL | w background | w pathogenic | direction vs manuscript | 2dlnL | LRT p | dS tree length (two-ratio) | max branch t |")
say("|---|---|---|---|---|---|---|---|---|---|---|")

summary = []
for key, (cfnum, cfname) in CF.items():
    m0 = runs.get(f"m0_{key}", (None, None))[0]
    tr = runs.get(f"two_ratio_{key}", (None, None))[0]
    if not (m0 and tr):
        say(f"| {cfnum} ({cfname}) | " + " | ".join(["pending"] * 10) + " |")
        continue
    d = 2 * (tr["lnL"] - m0["lnL"])
    p = chi2.sf(d, 1)
    direction = ("**supports** (pathogenic lower)" if tr["w_fg"] < tr["w_bg"]
                 else "**CONTRADICTS** (pathogenic higher)")
    bound = " **(at codeml bound 50)**" if tr["max_t"] and tr["max_t"] >= 49.99 else ""
    say(f"| {cfnum} ({cfname}) | {m0['lnL']:.3f} | {m0['omega']:.5f} | "
        f"{tr['lnL']:.3f} | {tr['w_bg']:.5f} | {tr['w_fg']:.5f} | {direction} | "
        f"{d:.3f} | {p:.3g} | {tr['dS_len']:.2f} | {tr['max_t']:.3f}{bound} |")
    summary.append((cfnum, cfname, m0, tr, d, p))

say("")
say("## Pending runs")
for name, (d, path) in sorted(runs.items()):
    if d is None:
        say(f"- `{path.parent.name}/{name}` — still running, no lnL in out.txt yet")
say("")

say("## Model fit between codon-frequency models")
say("")
say("PAML does not count codon frequencies in `np`, so lnL alone flatters the "
    "richer frequency models. Counting them explicitly (equal = 0 free, F1x4 = 3, "
    "F3x4 = 9, F61 = 60):")
say("")
say("| CodonFreq | M0 lnL | extra freq params | AIC (M0) | dAIC vs best |")
say("|---|---|---|---|---|")
extra = {"0": 0, "1": 3, "2": 9, "3": 60}
aics = []
for cfnum, cfname, m0, tr, d, p in summary:
    aic = 2 * (m0["np"] + extra[cfnum]) + 2 * (-m0["lnL"])
    aics.append((cfnum, cfname, m0, aic))
best = min(a[3] for a in aics)
for cfnum, cfname, m0, aic in aics:
    say(f"| {cfnum} ({cfname}) | {m0['lnL']:.3f} | {extra[cfnum]} | {aic:.1f} | {aic-best:+.1f} |")
say("")

(OUT / "R4.3_codonfreq_sweep_summary.md").write_text("\n".join(rep) + "\n")
print("\nwrote", OUT / "R4.3_codonfreq_sweep_summary.md")
