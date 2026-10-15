#!/usr/bin/env python3
"""Compare the manuscript's pairwise estimator against model-based estimates.

Addresses R3.2 ("add a model-based estimate ... and show that the ordering of
the three comparison groups is unchanged") and R2.2.

Three estimators on the SAME 22 species:
  1. manuscript   pn/ps, uncorrected proportions   (Work1/Genus_Diversity)
  2. NG86         Nei-Gojobori with Jukes-Cantor correction   (yn00 section A)
  3. YN00         Yang & Nielsen 2000: corrects ts/tv bias AND codon usage
                  (yn00 section B) -- the estimator R3.2 asks for

Also reports saturation diagnostics, since R3.2's specific claim is that the
deepest (mostly environmental-environmental) comparisons are least reliable.
"""
import sys, re, csv, math, itertools
from pathlib import Path
from statistics import mean, median

RES = Path(sys.argv[1])          # results dir (contains yn00/)
INP = Path(sys.argv[2])          # inputs dir (labels.tsv)
GD = Path(sys.argv[3])           # Work1/Genus_Diversity
OUT = Path(sys.argv[4])

labels = {r["abbrev"]: r for r in csv.DictReader((INP / "labels.tsv").open(), delimiter="\t")}
order = list(labels)             # alignment order == phylip order == yn00 indexing
cat = {k: v["category"] for k, v in labels.items()}


def group_of(a, b):
    ca, cb = cat[a], cat[b]
    if ca == cb == "pathogenic":
        return "pathogenic-pathogenic"
    if ca == cb == "environmental":
        return "environmental-environmental"
    return "cross"


BAD = {-1.0, 99.0, 999.0, -99.0}


def clean(x):
    """PAML writes -1 / 99 / 999 for estimates it could not compute."""
    return None if (x is None or x in BAD or not math.isfinite(x)) else x


# ---------------------------------------------------------------- YN00 table
yn = {}
txt = (RES / "yn00" / "yn00.out").read_text()
sec_b = txt.split("(B) Yang & Nielsen (2000) method", 1)[1]
pat = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+"
    r"([-\d.]+)\s+([-\d.]+)\s*\+-\s*([-\d.]+)\s+([-\d.]+)\s*\+-\s*([-\d.]+)",
    re.M)
for m in pat.finditer(sec_b):
    i, j = int(m.group(1)) - 1, int(m.group(2)) - 1
    a, b = order[i], order[j]
    yn[frozenset((a, b))] = {
        "S": float(m.group(3)), "N": float(m.group(4)),
        "t": float(m.group(5)), "kappa": float(m.group(6)),
        "omega": clean(float(m.group(7))),
        "dN": clean(float(m.group(8))), "dS": clean(float(m.group(10))),
    }

# ---------------------------------------------------------------- NG86 table
ng = {}
sec_a = txt.split("(A) Nei-Gojobori (1986) method", 1)[1].split("(B) Yang & Nielsen", 1)[0]
cellpat = re.compile(r"([-\d.]+)\s*\(\s*([-\d.]+)\s+([-\d.]+)\s*\)")
rows = [ln for ln in sec_a.splitlines() if cellpat.search(ln)]
for ln in rows:
    name = ln.split()[0]
    if name not in labels:
        continue
    i = order.index(name)
    cells = cellpat.findall(ln)
    for j, (om, dn, ds) in enumerate(cells):
        ng[frozenset((order[i], order[j]))] = {
            "omega": clean(float(om)), "dN": clean(float(dn)), "dS": clean(float(ds))}

# ------------------------------------------------------- manuscript estimates
ms = {}
with (GD / "results" / "pairwise_dnds_matrix.tsv").open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        ms[frozenset((r["species_A"], r["species_B"]))] = {
            "omega": float(r["dNdS"]), "ps": float(r["ps"]), "pn": float(r["pn"])}

# ------------------------------------------------------------------ assemble
pairs = []
for a, b in itertools.combinations(order, 2):
    k = frozenset((a, b))
    pairs.append({
        "species_A": a, "species_B": b, "comparison_type": group_of(a, b),
        "ms_omega": ms.get(k, {}).get("omega"), "ms_ps": ms.get(k, {}).get("ps"),
        "ng_omega": ng.get(k, {}).get("omega"), "ng_dS": ng.get(k, {}).get("dS"),
        "yn_omega": yn.get(k, {}).get("omega"), "yn_dS": yn.get(k, {}).get("dS"),
        "yn_kappa": yn.get(k, {}).get("kappa"),
    })

OUT.mkdir(parents=True, exist_ok=True)
with (OUT / "pairwise_three_estimators.tsv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(pairs[0]))
    w.writeheader()
    w.writerows(pairs)

GROUPS = ["pathogenic-pathogenic", "cross", "environmental-environmental"]


def summarise(key):
    out = {}
    for g in GROUPS:
        vals = [p[key] for p in pairs if p["comparison_type"] == g and p[key] is not None]
        miss = sum(1 for p in pairs if p["comparison_type"] == g and p[key] is None)
        out[g] = {"n": len(vals), "missing": miss,
                  "mean": mean(vals) if vals else None,
                  "median": median(vals) if vals else None}
    return out


def kruskal(key):
    """Kruskal-Wallis H with tie correction; chi-square survival via series."""
    data = [[p[key] for p in pairs if p["comparison_type"] == g and p[key] is not None]
            for g in GROUPS]
    allv = [v for d in data for v in d]
    N = len(allv)
    if min(len(d) for d in data) < 2:
        return None, None, N
    ranks = {}
    srt = sorted(range(N), key=lambda i: allv[i])
    i = 0
    while i < N:
        j = i
        while j + 1 < N and allv[srt[j + 1]] == allv[srt[i]]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[srt[k]] = r
        i = j + 1
    idx, H = 0, 0.0
    ties = {}
    for v in allv:
        ties[v] = ties.get(v, 0) + 1
    for d in data:
        rs = sum(ranks[idx + k] for k in range(len(d)))
        idx += len(d)
        H += rs ** 2 / len(d)
    H = 12.0 / (N * (N + 1)) * H - 3 * (N + 1)
    T = sum(t ** 3 - t for t in ties.values())
    if T:
        H /= (1 - T / (N ** 3 - N))
    # p-value for df=2 is closed form: exp(-H/2)
    return H, math.exp(-H / 2), N


lines = []
A = lines.append
A("# R3.2 / R2.2 — model-based re-estimation of pairwise dN/dS")
A("")
A("Same 22 species, same codon alignment, three estimators.")
A("")
A("| estimator | what it corrects | source |")
A("|---|---|---|")
A("| manuscript `pn/ps` | **nothing** — raw uncorrected proportions | Work1/Genus_Diversity |")
A("| NG86 | multiple hits (Jukes-Cantor) | yn00 section A |")
A("| YN00 | multiple hits + **ts/tv bias** + **codon usage** | yn00 section B |")
A("")

for key, name in [("ms_omega", "manuscript pn/ps"), ("ng_omega", "NG86 (JC-corrected)"),
                  ("yn_omega", "YN00 (ts/tv + codon usage)")]:
    s = summarise(key)
    H, p, N = kruskal(key)
    A(f"## {name}")
    A("")
    A("| comparison group | n | missing | mean | median |")
    A("|---|---|---|---|---|")
    for g in GROUPS:
        v = s[g]
        fmt = lambda x: "n/a" if x is None else f"{x:.4f}"
        A(f"| {g} | {v['n']} | {v['missing']} | {fmt(v['mean'])} | {fmt(v['median'])} |")
    A("")
    present = [(g, s[g]["mean"]) for g in GROUPS if s[g]["mean"] is not None]
    rank = sorted(present, key=lambda t: t[1])
    A("Ordering by mean: " + " < ".join(f"**{g}** ({v:.4f})" for g, v in rank))
    if H is not None:
        A(f"Kruskal-Wallis H = {H:.3f}, df = 2, p = {p:.3e} (n = {N})")
    A("")

# -------------------------------------------------------- saturation section
A("## Saturation diagnostics (R3.2's specific concern)")
A("")
A("| group | mean manuscript ps | mean NG86 dS | mean YN00 dS | NG86 undefined | YN00 undefined |")
A("|---|---|---|---|---|---|")
for g in GROUPS:
    sub = [p for p in pairs if p["comparison_type"] == g]
    def m(k):
        v = [p[k] for p in sub if p[k] is not None]
        return f"{mean(v):.4f}" if v else "n/a"
    A(f"| {g} | {m('ms_ps')} | {m('ng_dS')} | {m('yn_dS')} | "
      f"{sum(1 for p in sub if p['ng_omega'] is None)}/{len(sub)} | "
      f"{sum(1 for p in sub if p['yn_omega'] is None)}/{len(sub)} |")
A("")
kap = [p["yn_kappa"] for p in pairs if p["yn_kappa"] is not None]
A(f"Estimated kappa (ts/tv) across pairs: mean {mean(kap):.3f}, "
  f"median {median(kap):.3f}, range {min(kap):.3f}-{max(kap):.3f}")
A("")
A("kappa far from 1 means the manuscript's estimator, which assumes no ts/tv")
A("bias, was mis-specified; R3.2 and R4.3 both make this point.")
A("")

# ps above the Jukes-Cantor limit
over = [p for p in pairs if p["ms_ps"] is not None and p["ms_ps"] >= 0.75]
A(f"## Pairs at or beyond the Jukes-Cantor limit (ps >= 0.75): {len(over)}/{len(pairs)}")
A("")
A("At ps >= 0.75 the Jukes-Cantor correction -3/4 ln(1 - 4/3 ps) is undefined:")
A("synonymous divergence is saturated and dS cannot be recovered by that")
A("correction at all. The manuscript's uncorrected pn/ps sidesteps this only by")
A("not attempting the correction.")
A("")
for g in GROUPS:
    sub = [p for p in pairs if p["comparison_type"] == g]
    o = [p for p in sub if p["ms_ps"] is not None and p["ms_ps"] >= 0.75]
    A(f"- {g}: {len(o)}/{len(sub)} pairs")
A("")

(OUT / "R3.2_estimator_comparison.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
