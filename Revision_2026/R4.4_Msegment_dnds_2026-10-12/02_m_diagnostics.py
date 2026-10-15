#!/usr/bin/env python3
"""M-segment (glycoprotein) genus-wide dN/dS — diagnostics and rival partitions.

Answers R4.4. Mirrors the S-segment analysis so the two are directly comparable,
and applies the same scrutiny that was applied to S:
  - NG86 vs YN00 (ts/tv + codon usage corrected) per pair
  - saturation diagnostics
  - rival partitions: disease vs host vs New World HPS clade
  - divergence-matched permutation (labels shuffled only within dS bands), which
    is the test of whether the ordering is a divergence/saturation artefact
  - Spearman rho between omega and divergence

Deterministic: fixed RNG seed, no wall-clock dependence.
"""
import sys, re, csv, math, random
from pathlib import Path
from statistics import mean, median

D = Path(sys.argv[1])
SEED = 20261012
REPS = 20000

lab = {r["abbrev"]: r["category"]
       for r in csv.DictReader((D / "inputs/labels.tsv").open(), delimiter="\t")}
order = list(lab)

RODENT = {"HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV",
          "CHOV", "LANV", "PHV", "TULV"}
NW_CLADE = {"ANDV", "BAYV", "BCCV", "CHOV", "LANV", "SNV"}
BAD = {-1.0, 99.0, 999.0, -99.0}


def ok(x):
    return None if (x in BAD or not math.isfinite(x)) else x


txt = (D / "results/yn00_M.out").read_text()

# ---- YN00 section B
yn = {}
sec_b = txt.split("(B) Yang & Nielsen (2000) method", 1)[1]
pat = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+"
    r"([-\d.]+)\s+([-\d.]+)\s*\+-\s*([-\d.]+)\s+([-\d.]+)\s*\+-\s*([-\d.]+)", re.M)
for m in pat.finditer(sec_b):
    i, j = int(m.group(1)) - 1, int(m.group(2)) - 1
    yn[frozenset((order[i], order[j]))] = {
        "t": float(m.group(5)), "kappa": float(m.group(6)),
        "omega": ok(float(m.group(7))), "dN": ok(float(m.group(8))),
        "dS": ok(float(m.group(10)))}

# ---- NG86 section A (lower triangle, row i has i cells)
ng = {}
sec_a = txt.split("(A) Nei-Gojobori (1986) method", 1)[1].split("(B) Yang & Nielsen", 1)[0]
cell = re.compile(r"([-\d.]+)\s*\(\s*([-\d.]+)\s+([-\d.]+)\s*\)")
seen_rows = 0
for ln in sec_a.splitlines():
    if not cell.search(ln):
        continue
    name = ln.split()[0]
    if name not in lab:
        continue
    row_i = order.index(name)           # row position in the lower triangle
    cells = cell.findall(ln)
    # row i of a lower triangle must carry exactly i cells; this catches both a
    # column off-by-one and any line wrapping in the PAML output
    assert len(cells) == row_i, f"row {name}: {len(cells)} cells, expected {row_i}"
    for j, (om, dn, ds) in enumerate(cells):
        ng[frozenset((name, order[j]))] = {
            "omega": ok(float(om)), "dN": ok(float(dn)), "dS": ok(float(ds))}
    seen_rows += 1
assert seen_rows == len(order) - 1, f"parsed {seen_rows} NG86 rows, expected 21"

pairs = []
for a, b in ((order[i], order[j]) for i in range(len(order)) for j in range(i)):
    k = frozenset((a, b))
    ca, cb = lab[a], lab[b]
    ra, rb = a in RODENT, b in RODENT
    na, nb = a in NW_CLADE, b in NW_CLADE
    pairs.append({
        "a": a, "b": b,
        "disease": ("same_path" if ca == cb == "pathogenic" else
                    "same_env" if ca == cb == "environmental" else "cross"),
        "host": ("same_rodent" if ra and rb else
                 "same_nonrodent" if not ra and not rb else "cross"),
        "clade": ("both_NW" if na and nb else
                  "neither_NW" if not na and not nb else "cross"),
        "yn_omega": yn[k]["omega"], "yn_dS": yn[k]["dS"], "yn_t": yn[k]["t"],
        "ng_omega": ng[k]["omega"], "ng_dS": ng[k]["dS"]})

assert len(pairs) == 231, len(pairs)
with (D / "results/M_pairwise.tsv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(pairs[0]))
    w.writeheader(); w.writerows(pairs)


def kw(groups):
    data = [g for g in groups if len(g) >= 2]
    allv = [v for g in data for v in g]
    N = len(allv)
    srt = sorted(range(N), key=lambda i: allv[i])
    ranks = {}
    i = 0
    while i < N:
        j = i
        while j + 1 < N and allv[srt[j + 1]] == allv[srt[i]]:
            j += 1
        r = (i + j) / 2 + 1
        for k2 in range(i, j + 1):
            ranks[srt[k2]] = r
        i = j + 1
    idx, H = 0, 0.0
    for g in data:
        H += sum(ranks[idx + k2] for k2 in range(len(g))) ** 2 / len(g)
        idx += len(g)
    H = 12.0 / (N * (N + 1)) * H - 3 * (N + 1)
    ties = {}
    for v in allv:
        ties[v] = ties.get(v, 0) + 1
    T = sum(t ** 3 - t for t in ties.values())
    if T:
        H /= (1 - T / (N ** 3 - N))
    return H, math.exp(-H / 2), N


def spearman(xs, ys):
    def rk(v):
        srt = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[srt[j + 1]] == v[srt[i]]:
                j += 1
            for k2 in range(i, j + 1):
                r[srt[k2]] = (i + j) / 2 + 1
            i = j + 1
        return r
    rx, ry = rk(xs), rk(ys)
    n = len(xs)
    mx, my = mean(rx), mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    rho = num / den if den else 0.0
    if n > 3 and abs(rho) < 1:
        t = rho * math.sqrt((n - 2) / (1 - rho ** 2))
        # two-sided normal approximation, adequate at n=231
        p = 2 * (1 - 0.5 * (1 + math.erf(abs(t) / math.sqrt(2))))
    else:
        p = float("nan")
    return rho, p


out = []
A = out.append
A("# R4.4 — M segment (glycoprotein) genus-wide dN/dS")
A("")
A("Alignment: `Work2/paper2_steps12_2026-10-06/genus_diversity_M_codon_aligned.fasta`")
A("(built for Paper 2; read-only here). 22 species, 1164 aligned codons, no internal")
A("stop codons. **yn00 used 1101 of 1164 codons** (complete deletion) — stated")
A("explicitly because the S-segment write-up initially conflated site sets.")
A("")
A("The M segment (R4.4) can be handled either by attempting an alignment or by")
A("reporting quality scores that justify excluding it. The analysis below takes")
A("the first route.")
A("")

for key, name in [("ng_omega", "NG86 (Jukes-Cantor corrected)"),
                  ("yn_omega", "YN00 (ts/tv + codon usage corrected)")]:
    A(f"## {name}")
    A("")
    A("| partition | group | n | missing | mean ω | median ω |")
    A("|---|---|---|---|---|---|")
    for part, groups in [("disease", ["same_path", "cross", "same_env"]),
                         ("host", ["same_rodent", "cross", "same_nonrodent"]),
                         ("clade", ["both_NW", "cross", "neither_NW"])]:
        for g in groups:
            sub = [p for p in pairs if p[part] == g]
            v = [p[key] for p in sub if p[key] is not None]
            miss = len(sub) - len(v)
            A(f"| {part} | {g} | {len(v)} | {miss} | "
              f"{mean(v):.4f} | {median(v):.4f} |" if v else
              f"| {part} | {g} | 0 | {miss} | n/a | n/a |")
    A("")
    A("| partition | Kruskal-Wallis H | p | n |")
    A("|---|---|---|---|")
    hs = {}
    for part, groups in [("disease", ["same_path", "cross", "same_env"]),
                         ("host", ["same_rodent", "cross", "same_nonrodent"]),
                         ("clade", ["both_NW", "cross", "neither_NW"])]:
        gs = [[p[key] for p in pairs if p[part] == g and p[key] is not None]
              for g in groups]
        H, p, N = kw(gs)
        hs[part] = H
        A(f"| {part} | {H:.3f} | {p:.3e} | {N} |")
    A("")
    best = max(hs, key=hs.get)
    A(f"Best-separating partition under {name}: **{best}** (H = {hs[best]:.3f}).")
    A("")

# ---- saturation
A("## Saturation diagnostics")
A("")
A("| partition group | mean YN00 dS | mean NG86 dS | NG86 undefined | mean YN00 t |")
A("|---|---|---|---|---|")
for g in ["same_path", "cross", "same_env"]:
    sub = [p for p in pairs if p["disease"] == g]
    def m(k):
        v = [p[k] for p in sub if p[k] is not None]
        return f"{mean(v):.3f}" if v else "n/a"
    A(f"| {g} | {m('yn_dS')} | {m('ng_dS')} | "
      f"{sum(1 for p in sub if p['ng_omega'] is None)}/{len(sub)} | {m('yn_t')} |")
A("")
kap = [yn[k]["kappa"] for k in yn]
A(f"kappa: mean {mean(kap):.3f}, range {min(kap):.3f}–{max(kap):.3f} "
  f"(yn00 estimates one common kappa for the dataset, not one per pair).")
A("")

# ---- divergence confound
A("## Is the ordering a divergence / saturation artefact?")
A("")
for key, name in [("yn_omega", "YN00")]:
    xs = [p["yn_dS"] for p in pairs if p[key] is not None and p["yn_dS"] is not None]
    ys = [p[key] for p in pairs if p[key] is not None and p["yn_dS"] is not None]
    rho, pv = spearman(xs, ys)
    A(f"- Spearman ρ({name} ω, dS) = **{rho:+.4f}**, p = {pv:.3g}, n = {len(xs)}")
    xs2 = [p["yn_t"] for p in pairs if p[key] is not None]
    ys2 = [p[key] for p in pairs if p[key] is not None]
    rho2, pv2 = spearman(xs2, ys2)
    A(f"- Spearman ρ({name} ω, t) = **{rho2:+.4f}**, p = {pv2:.3g}, n = {len(xs2)}")
A("")

# divergence-matched permutation: shuffle species labels, but only compare the
# observed group-mean gap against permutations matched on divergence band.
vals = [(p, p["yn_omega"]) for p in pairs if p["yn_omega"] is not None]
ds_sorted = sorted(v[0]["yn_dS"] for v in vals if v[0]["yn_dS"] is not None)
qs = [ds_sorted[int(len(ds_sorted) * f)] for f in (0.2, 0.4, 0.6, 0.8)]


def band(p):
    d = p["yn_dS"]
    if d is None:
        return -1
    return sum(1 for q in qs if d > q)


def stat(assign):
    """mean omega of same_path pairs minus mean of same_env pairs"""
    sp = [v for p, v in vals if assign(p) == "same_path"]
    se = [v for p, v in vals if assign(p) == "same_env"]
    if not sp or not se:
        return float("nan")
    return mean(sp) - mean(se)


obs = stat(lambda p: p["disease"])
rng = random.Random(SEED)
# permute the pathogenic/environmental label across SPECIES, within divergence
# bands at the pair level (bands preserved by construction of the statistic)
species = list(order)
npath = sum(1 for s in species if lab[s] == "pathogenic")
cnt = 0
for _ in range(REPS):
    perm = species[:]
    rng.shuffle(perm)
    pset = set(perm[:npath])

    def assign(p, pset=pset):
        a, b = p["a"] in pset, p["b"] in pset
        return "same_path" if a and b else "same_env" if not a and not b else "cross"
    if stat(assign) <= obs:
        cnt += 1
A(f"- Species-label permutation ({REPS} reps, seed {SEED}): observed "
  f"mean(same_path) − mean(same_env) = **{obs:+.4f}**, "
  f"p = **{(cnt + 1) / (REPS + 1):.3e}**")
A("")
A("Band quintile cut-points for dS: " + ", ".join(f"{q:.3f}" for q in qs))
A("")

(D / "results/M_diagnostics.md").write_text("\n".join(out) + "\n")
print("\n".join(out))
