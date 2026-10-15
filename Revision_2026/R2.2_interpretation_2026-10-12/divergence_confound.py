#!/usr/bin/env python3
"""Is the pathogenic-vs-environmental pn/ps difference explained by divergence?

One concern (R3.2) is that the pairwise estimator is saturation-biased; another
(R3.1) is that the design cannot separate the effect from confounders. Both are testable
directly on the existing per-pair table: pathogenic-pathogenic pairs are the
LEAST divergent pairs in the matrix, and the manuscript's statistic depends on
divergence. Read-only on the R2.2 table; nothing is recomputed.
"""
import csv, json, math, random, statistics, sys
TSV, DEST = sys.argv[1:3]
random.seed(20261012)
rows = []
for r in csv.DictReader(open(TSV), delimiter="\t"):
    try:
        rows.append({"grp": r["comparison_type"], "om": float(r["ms_omega"]),
                     "ps": float(r["ms_ps"]), "yn_om": float(r["yn_omega"]),
                     "yn_dS": float(r["yn_dS"])})
    except ValueError:
        continue

def spearman(x, y):
    def rank(v):
        s = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
        while i < len(s):
            j = i
            while j+1 < len(s) and v[s[j+1]] == v[s[i]]: j += 1
            for k in range(i, j+1): r[s[k]] = (i+j)/2.0+1
            i = j+1
        return r
    rx, ry = rank(x), rank(y); n = len(x)
    mx, my = sum(rx)/n, sum(ry)/n
    num = sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    den = math.sqrt(sum((v-mx)**2 for v in rx)*sum((v-my)**2 for v in ry))
    rho = num/den
    t = rho*math.sqrt((n-2)/max(1e-12, 1-rho*rho))
    # two-sided p from normal approx on Fisher z
    z = 0.5*math.log((1+rho)/(1-rho))*math.sqrt(n-3)
    return rho, math.erfc(abs(z)/math.sqrt(2))

GRPS = ["pathogenic-pathogenic", "cross", "environmental-environmental"]
out = {"n_pairs": len(rows)}
out["divergence_by_group_ms_ps"] = {g: {
    "n": sum(1 for r in rows if r["grp"] == g),
    "mean_ps": statistics.mean([r["ps"] for r in rows if r["grp"] == g]),
    "mean_yn_dS": statistics.mean([r["yn_dS"] for r in rows if r["grp"] == g]),
} for g in GRPS}
rho_ps, p_ps = spearman([r["ps"] for r in rows], [r["om"] for r in rows])
rho_ds, p_ds = spearman([r["yn_dS"] for r in rows], [r["yn_om"] for r in rows])
out["spearman_manuscript_omega_vs_ps"] = {"rho": rho_ps, "p": p_ps}
out["spearman_yn_omega_vs_yn_dS"] = {"rho": rho_ds, "p": p_ds}

# stratify on divergence: within narrow ps bands, does the group ordering hold?
bands = [(0.0, 0.76), (0.76, 0.79), (0.79, 0.81), (0.81, 1.01)]
strat = []
for lo, hi in bands:
    sub = [r for r in rows if lo <= r["ps"] < hi]
    rec = {"ps_band": [lo, hi], "n": len(sub)}
    for g in GRPS:
        v = [r["om"] for r in sub if r["grp"] == g]
        rec[g] = {"n": len(v), "mean": statistics.mean(v) if v else None}
    ok = [rec[g]["mean"] for g in GRPS if rec[g]["mean"] is not None]
    rec["ordering_preserved"] = (len(ok) == 3 and ok[0] < ok[1] < ok[2])
    rec["testable"] = all(rec[g]["n"] >= 5 for g in GRPS)
    strat.append(rec)
out["stratified_by_divergence"] = strat

# divergence-matched permutation: permute group labels among pairs of SIMILAR
# divergence only (within ps band), so the null keeps the divergence structure
def group_means(labels):
    d = {g: [] for g in GRPS}
    for r, g in zip(rows, labels):
        d[g].append(r["om"])
    return {g: statistics.mean(v) if v else None for g, v in d.items()}
obs = group_means([r["grp"] for r in rows])
obs_gap = obs["environmental-environmental"] - obs["pathogenic-pathogenic"]
idx_by_band = {}
for i, r in enumerate(rows):
    for lo, hi in bands:
        if lo <= r["ps"] < hi:
            idx_by_band.setdefault((lo, hi), []).append(i)
            break
REPS = 20000
ge = 0
for _ in range(REPS):
    lab = [None]*len(rows)
    for band, idxs in idx_by_band.items():
        labs = [rows[i]["grp"] for i in idxs]
        random.shuffle(labs)
        for i, g in zip(idxs, labs):
            lab[i] = g
    m = group_means(lab)
    gap = m["environmental-environmental"] - m["pathogenic-pathogenic"]
    if gap >= obs_gap:
        ge += 1
out["divergence_matched_permutation"] = {
    "reps": REPS, "observed_env_minus_path": obs_gap,
    "p_one_sided": (ge+1)/(REPS+1),
    "note": "labels permuted only WITHIN divergence bands, so the null retains "
            "the divergence structure that the naive permutation test destroys",
}
json.dump(out, open(DEST+"/divergence_confound.json", "w"), indent=2, sort_keys=True)
print(json.dumps(out, indent=2))
