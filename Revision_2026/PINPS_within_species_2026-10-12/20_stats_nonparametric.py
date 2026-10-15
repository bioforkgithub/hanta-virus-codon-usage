#!/usr/bin/env python3
"""
Nonparametric / exact statistical treatment of the per-species within-species
piN/piS table produced by phase 2 of PINPS_within_species_2026-10-12.

Everything here is exact (complete enumeration) where exact is possible, and
simulation with a fixed seed where it is not.  Nothing is estimated by eye and
nothing is carried over from memory: the only input is
    <results>/pinps_per_species.tsv
plus the species labels in
    R2.2_phylo_dnds_2026-10-12/inputs/labels.tsv
and, for the single-clade relabelling, the ML tree
    R2.2_phylo_dnds_2026-10-12/inputs/S22_ml_clean.nwk

Usage:
    python3 -I 20_stats_nonparametric.py <PINPS_ANALYSIS_ROOT> <R22_ANALYSIS_ROOT> <OUTDIR>

Deterministic.  Master seed MASTER_SEED below; no wall-clock input anywhere.
"""

import sys
import os
import json
import math
import itertools
from collections import OrderedDict

MASTER_SEED = 20261012

# ---------------------------------------------------------------- utilities


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            parts += [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts)))
    return header, rows


def fnum(s):
    s = (s or "").strip()
    if s == "" or s.lower() in ("na", "nan", "none"):
        return None
    return float(s)


def median(xs):
    ys = sorted(xs)
    n = len(ys)
    if n == 0:
        return None
    if n % 2:
        return ys[n // 2]
    return 0.5 * (ys[n // 2 - 1] + ys[n // 2])


def mean(xs):
    return sum(xs) / len(xs) if xs else None


def quantile(xs, q):
    """Type-7 (R default) quantile."""
    ys = sorted(xs)
    n = len(ys)
    if n == 0:
        return None
    if n == 1:
        return ys[0]
    h = (n - 1) * q
    lo = int(math.floor(h))
    hi = min(lo + 1, n - 1)
    return ys[lo] + (h - lo) * (ys[hi] - ys[lo])


# ------------------------------------------------- exact Mann-Whitney / W

def midranks(values):
    """Midranks (average ranks for ties), 1-based."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        r = 0.5 * ((i + 1) + (j + 1))
        for k in range(i, j + 1):
            ranks[order[k]] = r
        i = j + 1
    return ranks


def exact_mwu(x, y):
    """
    Exact two-sample rank-sum test by complete enumeration of the
    C(n, nx) label assignments of the OBSERVED (mid)ranks.  This is the exact
    conditional permutation distribution of the rank-sum statistic and is
    correct in the presence of ties.

    x = group 1 values, y = group 2 values.
    Returns a dict.  U1 counts (number of (xi, yj) pairs with xi > yj) plus
    half the pairs that tie, i.e. the standard Mann-Whitney U for group 1.
    """
    nx, ny = len(x), len(y)
    n = nx + ny
    allv = list(x) + list(y)
    R = midranks(allv)
    Rx = sum(R[:nx])
    U1 = Rx - nx * (nx + 1) / 2.0          # "x beats y" count (+0.5 per tie)
    U2 = nx * ny - U1
    # exact null distribution of Rx
    idx = list(range(n))
    tot = 0
    ge = 0
    le = 0
    extreme = 0
    dist = {}
    for comb in itertools.combinations(idx, nx):
        s = sum(R[i] for i in comb)
        dist[s] = dist.get(s, 0) + 1
        tot += 1
    mu = nx * (n + 1) / 2.0
    for s, c in dist.items():
        if s >= Rx - 1e-12:
            ge += c
        if s <= Rx + 1e-12:
            le += c
        if abs(s - mu) >= abs(Rx - mu) - 1e-12:
            extreme += c
    p_two = extreme / tot
    # one-sided: "x stochastically smaller than y"  <=>  small Rx
    p_x_lt_y = le / tot
    p_x_gt_y = ge / tot
    # minimum attainable p-values for this design (and these ties)
    smin = min(dist)
    smax = max(dist)
    nmin = dist[smin]
    nmax = dist[smax]
    p_two_min = (nmin + nmax) / tot if abs(smin - mu) == abs(smax - mu) else \
        max(nmin, nmax) / tot
    p_one_min = min(nmin, nmax) / tot
    return {
        "n_x": nx, "n_y": ny,
        "rank_sum_x_W": Rx,
        "rank_sum_x_null_mean": mu,
        "U_x_over_y": U1,
        "U_y_over_x": U2,
        "theta_hat_P(y>x)": U2 / (nx * ny) if nx * ny else None,
        "n_permutations_enumerated": tot,
        "p_two_sided_exact": p_two,
        "p_one_sided_exact_x_lt_y": p_x_lt_y,
        "p_one_sided_exact_x_gt_y": p_x_gt_y,
        "min_attainable_p_two_sided": p_two_min,
        "min_attainable_p_one_sided": p_one_min,
        "has_ties": len(set(allv)) != len(allv),
    }


def hodges_lehmann(x, y):
    """HL estimate of the shift y - x (median of all pairwise differences),
    and the exact distribution-free CI from the exact null of U."""
    diffs = sorted(yy - xx for xx in x for yy in y)
    hl = median(diffs)
    nx, ny = len(x), len(y)
    m = nx * ny
    # exact null distribution of U (no ties assumed for the CI; U integer)
    rs = ranksum_null(nx + ny, nx)
    off = nx * (nx + 1) // 2
    dist = {s - off: c for s, c in rs.items()}
    tot = sum(dist.values())
    # largest k such that P(U <= k-1) <= 0.025
    cum = 0
    k = 0
    for u in range(0, m + 1):
        c = dist.get(u, 0)
        if (cum + c) / tot > 0.025:
            break
        cum += c
        k = u + 1
    if k == 0:
        return {"HL_shift_y_minus_x": hl, "ci_level_achieved": None,
                "ci_lo": None, "ci_hi": None,
                "ci_note": "no distribution-free CI of level >=95% exists at this n"}
    lo = diffs[k - 1]
    hi = diffs[m - k]
    cover = 1 - 2 * cum / tot
    return {"HL_shift_y_minus_x": hl,
            "ci_level_achieved": cover,
            "ci_lo": lo, "ci_hi": hi,
            "ci_k_order_statistic": k,
            "ci_note": "exact distribution-free (Moses) interval, no-tie null"}


def exact_permutation(x, y, stats=("diff_median", "diff_mean", "ratio_mean")):
    """Exact permutation test of the group labels by complete enumeration of
    all C(n, nx) assignments of the observed VALUES (not ranks)."""
    nx, ny = len(x), len(y)
    allv = list(x) + list(y)
    n = nx + ny
    idx = list(range(n))

    def compute(xs, ys):
        out = {}
        out["diff_median"] = median(ys) - median(xs)
        out["diff_mean"] = mean(ys) - mean(xs)
        out["ratio_mean"] = (mean(ys) / mean(xs)) if mean(xs) > 0 else float("inf")
        return out

    obs = compute(x, y)
    counts = {s: 0 for s in stats}
    tot = 0
    nulls = {s: [] for s in stats}
    for comb in itertools.combinations(idx, nx):
        cs = set(comb)
        xs = [allv[i] for i in comb]
        ys = [allv[i] for i in idx if i not in cs]
        v = compute(xs, ys)
        tot += 1
        for s in stats:
            nulls[s].append(v[s])
    res = {"n_permutations_enumerated": tot}
    for s in stats:
        o = obs[s]
        centre = 1.0 if s == "ratio_mean" else 0.0
        if s == "ratio_mean":
            # symmetrise on the log scale
            dev = abs(math.log(o)) if o > 0 and math.isfinite(o) else float("inf")
            nd = [abs(math.log(v)) if v > 0 and math.isfinite(v) else float("inf")
                  for v in nulls[s]]
        else:
            dev = abs(o - centre)
            nd = [abs(v - centre) for v in nulls[s]]
        ge = sum(1 for v in nd if v >= dev - 1e-15)
        res[s] = {
            "observed": o,
            "p_two_sided_exact": ge / tot,
            "p_one_sided_exact_y_gt_x": sum(1 for v in nulls[s] if v >= o - 1e-15) / tot,
            "null_min": min(nulls[s]), "null_max": max(nulls[s]),
        }
    return res


# ----------------------------------------------------- power / MDE by sim

_RSN_CACHE = {}


def ranksum_null(n, k):
    """Exact counts of k-subsets of {1..n} by their sum (dynamic programming).
    Equivalent to full enumeration but feasible for large n."""
    key = (n, k)
    if key in _RSN_CACHE:
        return _RSN_CACHE[key]
    # dp[j] = dict sum -> count using j chosen elements
    dp = [dict() for _ in range(k + 1)]
    dp[0][0] = 1
    for v in range(1, n + 1):
        for j in range(min(k, v), 0, -1):
            src = dp[j - 1]
            dst = dp[j]
            for s, c in src.items():
                dst[s + v] = dst.get(s + v, 0) + c
    _RSN_CACHE[key] = dp[k]
    return dp[k]


def power_lehmann(k_path, k_env, lam, alpha, n_sim, rng, two_sided=True):
    """
    Distribution-free power of the exact rank-sum test under a Lehmann
    alternative: F_env(t) = F_path(t)**lam, realised as
        X_path ~ Uniform(0,1)             (so F_path = identity on [0,1])
        X_env  ~ Uniform(0,1)**(1/lam)    (so F_env(t) = t**lam)
    lam > 1  =>  environmental values stochastically LARGER
    lam < 1  =>  environmental values stochastically SMALLER
    theta = P(X_env > X_path) = integral_0^1 (1 - t**lam) dt = lam/(1+lam).
    Because the test is rank-based, the power depends on lam alone and not on
    the marginal distribution, so this is an assumption-free power curve in
    theta.  Rejection region is read off the exact null of the rank sum.
    """
    import numpy as np
    n = k_path + k_env
    # exact null of rank sum of the pathogenic group (untied ranks), by DP
    dist = ranksum_null(n, k_path)
    tot = sum(dist.values())
    mu = k_path * (n + 1) / 2.0
    # critical set: all rank sums whose two-sided (or one-sided) exact p <= alpha
    def p_two(s):
        return sum(c for t, c in dist.items()
                   if abs(t - mu) >= abs(s - mu) - 1e-12) / tot

    def p_one(s):
        return sum(c for t, c in dist.items() if t <= s + 1e-12) / tot

    crit = set()
    for s in dist:
        p = p_two(s) if two_sided else p_one(s)
        if p <= alpha + 1e-15:
            crit.add(s)
    if not crit:
        return {"power": 0.0, "theta": lam / (1.0 + lam),
                "critical_rank_sums": [], "note": "no rejection possible at this alpha"}
    u_p = rng.random((n_sim, k_path))
    u_e = rng.random((n_sim, k_env)) ** (1.0 / lam)
    allv = np.concatenate([u_p, u_e], axis=1)
    # rank sum of the first k_path columns
    order = np.argsort(allv, axis=1, kind="stable")
    ranks = np.empty_like(order)
    rows = np.arange(n_sim)[:, None]
    ranks[rows, order] = np.arange(1, n + 1)[None, :]
    rs = ranks[:, :k_path].sum(axis=1)
    hit = np.isin(rs, np.array(sorted(crit)))
    return {"power": float(hit.mean()), "theta": lam / (1.0 + lam),
            "lambda": lam,
            "critical_rank_sums": sorted(crit),
            "n_sim": n_sim}


def theta_to_fold(theta, sigma_log):
    """Translate a probability-of-superiority theta into a multiplicative
    fold-difference in the median, under an explicitly stated lognormal model
    with equal log-SD sigma_log in both groups:
        theta = Phi(Delta / (sigma_log * sqrt(2))),  fold = exp(Delta).
    """
    if not (0 < theta < 1):
        return None
    # inverse normal CDF (Acklam / Moro-free: use erfinv via math)
    z = math.sqrt(2.0) * _erfinv(2.0 * theta - 1.0)
    return math.exp(z * sigma_log * math.sqrt(2.0))


def _erfinv(y):
    """Inverse error function, Newton-refined from a rational start."""
    if y <= -1.0 or y >= 1.0:
        return float("-inf") if y <= -1 else float("inf")
    a = 0.147
    ln1 = math.log(1 - y * y)
    t1 = 2.0 / (math.pi * a) + ln1 / 2.0
    x = math.copysign(math.sqrt(max(0.0, math.sqrt(t1 * t1 - ln1 / a) - t1)), y)
    for _ in range(60):
        err = math.erf(x) - y
        d = 2.0 / math.sqrt(math.pi) * math.exp(-x * x)
        if d == 0:
            break
        x -= err / d
    return x


# ------------------------------------------------------------------- main

def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    pinps_root, r22_root, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)

    import numpy as np

    tsv = os.path.join(pinps_root, "results", "pinps_per_species.tsv")
    labels_p = os.path.join(r22_root, "inputs", "labels.tsv")
    tree_p = os.path.join(r22_root, "inputs", "S22_ml_clean.nwk")

    _, rows = read_tsv(tsv)
    _, lab = read_tsv(labels_p)
    labels = {r["abbrev"]: r for r in lab if r.get("abbrev")}
    with open(tree_p) as fh:
        newick = fh.read().strip()

    # ---- New World HPS clade membership, read off the ML tree, not assumed.
    # The clade is the smallest bipartition of S22_ml_clean.nwk containing SNV.
    nw_clade = nw_hps_clade(newick)

    est_cols = OrderedDict([
        ("piN_piS", "primary: (mean pN)/(mean pS) over pairs"),
        ("piN_piS_tot", "ratio of totals (additive over codons)"),
        ("piN_piS_jc", "Jukes-Cantor corrected"),
        ("piN_piS_stopexcl", "MEGA stop-excluding site convention"),
    ])

    # ---------- build the analysis panels
    panels = OrderedDict()

    def panel(name, seqset, level, gate, note):
        sel = []
        for r in rows:
            if r["seqset"] != seqset or r["n_level"] != level:
                continue
            if gate == "usable" and r["usable_for_group_test"] != "True":
                continue
            if gate == "n10":
                if fnum(r["n_sequences"]) is None or fnum(r["n_sequences"]) < 10:
                    continue
            sel.append(r)
        if sel:
            panels[name] = {"seqset": seqset, "n_level": level, "gate": gate,
                            "note": note, "rows": sel}

    panel("P1_clust95d_full", "clust95d", "full", "usable",
          "PRIMARY descriptive: single-population set, saturation+pool gates passed, "
          "each species at its own full n")
    panel("P2_clust95d_n10", "clust95d", "10", "usable",
          "PRIMARY effort-controlled: every species subsampled to n=10, 1000 draws. "
          "This is the headline panel.")
    panel("P3_clust95d_n12", "clust95d", "12", "usable",
          "equal-n ladder, n=12 (CHOV has only 12 sequences so it has no n=12 "
          "subsample row and drops out)")
    panel("P4_clust95d_n15", "clust95d", "15", "usable",
          "equal-n ladder, n=15")
    panel("P5_clust95_full", "clust95", "full", "usable",
          "clustering sensitivity: single-linkage cluster instead of the "
          "diameter-pruned one")
    panel("P6_ss_all_full_n10", "ss_all", "full", "n10",
          "PREMISE-VIOLATING sensitivity: whole ICTV species, n>=10, saturation gate "
          "IGNORED (11 of these 12 species have piS_jc ~ 1). Reported because it is "
          "the comparison a reader would naively make; it is not valid.")
    panel("P7_ss_all_n10", "ss_all", "10", "n10",
          "PREMISE-VIOLATING sensitivity at equal n=10, saturation gate IGNORED")

    out = {
        "meta": {
            "master_seed": MASTER_SEED,
            "input_table": tsv,
            "input_labels": labels_p,
            "input_tree": tree_p,
            "new_world_hps_clade_from_tree": sorted(nw_clade),
            "estimators_tested": list(est_cols.keys()),
            "note": "All p-values are EXACT (complete enumeration). No normal "
                    "approximation and no continuity correction is used anywhere.",
        },
        "panels": OrderedDict(),
    }

    labellings = OrderedDict([
        ("disease", lambda r: r["category"]),
        ("disease_TULV_moved", lambda r: r["category_tulv_moved"]),
        ("host_class_rodent", lambda r: "rodent" if labels[r["abbrev"]]["is_rodent"] == "True" else "nonrodent"),
        ("NW_HPS_clade", lambda r: "NWclade" if r["abbrev"] in nw_clade else "other"),
    ])
    # group1 is the "focal / expected-lower" group in each labelling
    focal = {"disease": "pathogenic", "disease_TULV_moved": "pathogenic",
             "host_class_rodent": "rodent", "NW_HPS_clade": "NWclade"}

    for pname, p in panels.items():
        rws = p["rows"]
        entry = {
            "definition": {k: p[k] for k in ("seqset", "n_level", "gate", "note")},
            "species": [],
        }
        for r in sorted(rws, key=lambda r: fnum(r["piN_piS"])):
            entry["species"].append({
                "abbrev": r["abbrev"],
                "category": r["category"],
                "category_tulv_moved": r["category_tulv_moved"],
                "host_class": r["host_class"],
                "is_rodent": labels[r["abbrev"]]["is_rodent"] == "True",
                "in_NW_HPS_clade": r["abbrev"] in nw_clade,
                "n_sequences": int(float(r["n_sequences"])),
                "n_countries": r["n_countries"],
                "n_human": r["n_human"],
                "flag_piS_saturated": r["flag_piS_saturated"],
                "piS_jc": fnum(r["piS_jc"]),
                "piN": fnum(r["piN"]),
                "piS": fnum(r["piS"]),
                "piN_piS": fnum(r["piN_piS"]),
                "piN_piS_tot": fnum(r["piN_piS_tot"]),
                "piN_piS_jc": fnum(r["piN_piS_jc"]),
                "piN_piS_stopexcl": fnum(r["piN_piS_stopexcl"]),
                "ratio_lo": fnum(r["ratio_lo"]),
                "ratio_hi": fnum(r["ratio_hi"]),
            })

        entry["tests"] = OrderedDict()
        for lname, lfun in labellings.items():
            g1 = [r for r in rws if lfun(r) == focal[lname]]
            g2 = [r for r in rws if lfun(r) != focal[lname]]
            key = "%s (%s=%d vs other=%d)" % (lname, focal[lname], len(g1), len(g2))
            sub = OrderedDict()
            sub["group1_label"] = focal[lname]
            sub["group1_species"] = [r["abbrev"] for r in g1]
            sub["group2_species"] = [r["abbrev"] for r in g2]
            if len(g1) == 0 or len(g2) == 0:
                sub["result"] = "UNDEFINED: one group is empty at this panel/labelling"
                entry["tests"][key] = sub
                continue
            for var in list(est_cols) + ["piS", "piN"]:
                x = [fnum(r[var]) for r in g1]
                y = [fnum(r[var]) for r in g2]
                if any(v is None for v in x + y):
                    continue
                t = OrderedDict()
                t["group1"] = {"n": len(x), "median": median(x), "mean": mean(x),
                               "min": min(x), "max": max(x),
                               "q25": quantile(x, .25), "q75": quantile(x, .75)}
                t["group2"] = {"n": len(y), "median": median(y), "mean": mean(y),
                               "min": min(y), "max": max(y),
                               "q25": quantile(y, .25), "q75": quantile(y, .75)}
                t["direction_median"] = ("group1 lower" if median(x) < median(y)
                                         else "group1 higher" if median(x) > median(y)
                                         else "equal")
                t["direction_mean"] = ("group1 lower" if mean(x) < mean(y)
                                       else "group1 higher" if mean(x) > mean(y)
                                       else "equal")
                t["mannwhitney_exact"] = exact_mwu(x, y)
                t["hodges_lehmann"] = hodges_lehmann(x, y)
                t["label_permutation_exact"] = exact_permutation(x, y)
                # CI overlap bookkeeping for the primary estimator only
                if var == "piN_piS":
                    ov = 0
                    for r1 in g1:
                        for r2 in g2:
                            a1, b1 = fnum(r1["ratio_lo"]), fnum(r1["ratio_hi"])
                            a2, b2 = fnum(r2["ratio_lo"]), fnum(r2["ratio_hi"])
                            if None in (a1, b1, a2, b2):
                                ov = None
                                break
                            if a1 <= b2 and a2 <= b1:
                                ov += 1
                        if ov is None:
                            break
                    t["cross_group_pairs_with_overlapping_95pct_intervals"] = ov
                    t["cross_group_pairs_total"] = len(g1) * len(g2)
                sub[var] = t
            entry["tests"][key] = sub

        # ---- restricted / deconfounded cells
        entry["restricted_cells"] = OrderedDict()
        for cname, filt in [
            ("rodent_borne_only__disease_labels",
             lambda r: labels[r["abbrev"]]["is_rodent"] == "True"),
            ("non_rodent_only__disease_labels",
             lambda r: labels[r["abbrev"]]["is_rodent"] != "True"),
            ("rodent_borne_only__TULV_moved",
             lambda r: labels[r["abbrev"]]["is_rodent"] == "True"),
        ]:
            lab_fn = (lambda r: r["category_tulv_moved"]) if cname.endswith("TULV_moved") \
                else (lambda r: r["category"])
            sel = [r for r in rws if filt(r)]
            g1 = [r for r in sel if lab_fn(r) == "pathogenic"]
            g2 = [r for r in sel if lab_fn(r) == "environmental"]
            c = {"n_pathogenic": len(g1), "n_environmental": len(g2),
                 "pathogenic_species": [r["abbrev"] for r in g1],
                 "environmental_species": [r["abbrev"] for r in g2]}
            if len(g1) == 0 or len(g2) == 0:
                c["result"] = ("UNDEFINED / NO TEST POSSIBLE: one group is empty "
                               "(this cell does not exist in the data)")
            else:
                x = [fnum(r["piN_piS"]) for r in g1]
                y = [fnum(r["piN_piS"]) for r in g2]
                c["pathogenic_values"] = x
                c["environmental_values"] = y
                c["mannwhitney_exact"] = exact_mwu(x, y)
                c["verdict"] = (
                    "UNDERPOWERED BY CONSTRUCTION: the smallest two-sided exact "
                    "p attainable with these group sizes is %.4f"
                    % c["mannwhitney_exact"]["min_attainable_p_two_sided"])
            entry["restricted_cells"][cname] = c

        # ---- rank invariance across estimators
        ranks_by_est = {}
        for var in est_cols:
            vals = [(r["abbrev"], fnum(r[var])) for r in rws]
            if any(v is None for _, v in vals):
                continue   # estimator not reported on subsample rows
            ranks_by_est[var] = [a for a, v in sorted(vals, key=lambda t: t[1])]
        entry["estimator_rank_order"] = ranks_by_est
        entry["estimator_rank_order_identical"] = (
            len(set(tuple(v) for v in ranks_by_est.values())) == 1)

        # ---- leave-one-species-out stability of the direction (disease labels)
        loo = []
        base = [r for r in rws]
        for drop in base:
            keep = [r for r in base if r is not drop]
            g1 = [fnum(r["piN_piS"]) for r in keep if r["category"] == "pathogenic"]
            g2 = [fnum(r["piN_piS"]) for r in keep if r["category"] == "environmental"]
            if not g1 or not g2:
                loo.append({"dropped": drop["abbrev"], "result": "a group becomes empty"})
                continue
            mw = exact_mwu(g1, g2)
            loo.append({
                "dropped": drop["abbrev"],
                "n_path": len(g1), "n_env": len(g2),
                "median_path": median(g1), "median_env": median(g2),
                "mean_path": mean(g1), "mean_env": mean(g2),
                "direction_median": "path lower" if median(g1) < median(g2) else "path higher",
                "direction_mean": "path lower" if mean(g1) < mean(g2) else "path higher",
                "p_two_sided_exact": mw["p_two_sided_exact"],
                "min_attainable_p_two_sided": mw["min_attainable_p_two_sided"],
            })
        entry["leave_one_species_out"] = loo

        # ---- exact permutation of the disease labels STRATIFIED by host class
        # (the only way to test disease while holding the host-class confound
        # fixed without discarding species).
        entry["host_class_stratified_permutation"] = stratified_perm(
            rws, labels, strat_key=lambda r: labels[r["abbrev"]]["host_class"])
        entry["rodent_vs_nonrodent_stratified_permutation"] = stratified_perm(
            rws, labels,
            strat_key=lambda r: ("rodent" if labels[r["abbrev"]]["is_rodent"] == "True"
                                 else "nonrodent"))

        # ---- the sister-pair (phylogenetically paired) contrast(s)
        entry["sister_pair_contrasts"] = sister_pairs(rws, labels, newick)

        # ---- nonparametric correlates: is the per-species value explained by
        # sampling depth, by clinical origin, or by geographic breadth?
        entry["spearman_correlates"] = spearman_block(rws, panels, pname)

        out["panels"][pname] = entry

    # ------------------------------------------------ power / MDE section
    rng = np.random.default_rng([MASTER_SEED, 99])
    N_SIM = 100000
    power = {"method": (
        "Distribution-free power of the EXACT rank-sum test under a Lehmann "
        "alternative F_env = F_path^lambda, so theta = P(env > path) = 1/(1+lambda). "
        "Because the test uses only ranks, this power curve depends on theta alone "
        "and on no distributional assumption. %d simulations per point, "
        "rng = numpy PCG64 seeded [%d, 99]."
        % (N_SIM, MASTER_SEED)),
        "designs": OrderedDict()}

    thetas = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.825, 0.85, 0.875,
              0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.955, 0.96, 0.965, 0.97,
              0.975, 0.98, 0.985, 0.99, 0.9925, 0.995, 0.9975, 0.999,
              0.9995, 0.9999]
    designs = [(7, 3), (6, 3), (8, 2), (7, 5), (5, 5), (10, 10), (12, 12),
               (15, 15), (20, 20), (25, 25), (30, 30)]
    for (kp, ke) in designs:
        d = {"k_path": kp, "k_env": ke, "curve_two_sided_alpha0.05": [],
             "curve_one_sided_alpha0.05": []}
        for th in thetas:
            lam = th / (1.0 - th)
            r2 = power_lehmann(kp, ke, lam, 0.05, N_SIM, rng, two_sided=True)
            r1 = power_lehmann(kp, ke, lam, 0.05, N_SIM, rng, two_sided=False)
            d["curve_two_sided_alpha0.05"].append({"theta": th, "power": r2["power"]})
            d["curve_one_sided_alpha0.05"].append({"theta": th, "power": r1["power"]})
        # theta needed for 80% power, by interpolation on the simulated curve
        for key in ("curve_two_sided_alpha0.05", "curve_one_sided_alpha0.05"):
            cur = d[key]
            need = None
            for a, b in zip(cur, cur[1:]):
                if a["power"] < 0.80 <= b["power"]:
                    f = (0.80 - a["power"]) / (b["power"] - a["power"])
                    need = a["theta"] + f * (b["theta"] - a["theta"])
                    break
            d["theta_for_80pct_power__" + key] = need
            if need is None:
                d["theta_for_80pct_power__" + key + "_note"] = (
                    "80%% power is NOT reached at any theta up to %.3f "
                    "(max simulated power %.4f)" % (thetas[-1], cur[-1]["power"]))
        power["designs"]["k%d_vs_k%d" % (kp, ke)] = d
    out["power_and_MDE"] = power

    # observed theta + fold translation for the headline panel
    hp = out["panels"].get("P2_clust95d_n10")
    if hp:
        k = [kk for kk in hp["tests"] if kk.startswith("disease ")][0]
        t = hp["tests"][k]["piN_piS"]
        xs = [s["piN_piS"] for s in hp["species"] if s["category"] == "pathogenic"]
        ys = [s["piN_piS"] for s in hp["species"] if s["category"] == "environmental"]
        pos = [v for v in xs if v > 0]
        sl = math.sqrt(sum((math.log(v) - mean([math.log(w) for w in pos])) ** 2
                           for v in pos) / (len(pos) - 1))
        mde = {
            "panel": "P2_clust95d_n10",
            "observed_theta_P(env>path)": t["mannwhitney_exact"]["theta_hat_P(y>x)"],
            "observed_p_two_sided_exact": t["mannwhitney_exact"]["p_two_sided_exact"],
            "min_attainable_p_two_sided": t["mannwhitney_exact"]["min_attainable_p_two_sided"],
            "min_attainable_p_one_sided": t["mannwhitney_exact"]["min_attainable_p_one_sided"],
            "observed_ratio_of_medians_env_over_path": median(ys) / median(xs),
            "observed_ratio_of_means_env_over_path": mean(ys) / mean(xs),
            "sigma_log_within_pathogenic_group": sl,
            "sigma_log_note": ("log-SD of the %d pathogenic piN/piS values with "
                               "piN/piS > 0 (SNV is exactly 0 and is excluded from "
                               "this one quantity only). Used ONLY to translate "
                               "theta into a fold-change for readability; the "
                               "p-values and the theta-scale MDE do not use it."
                               % len(pos)),
        }
        d73 = power["designs"]["k7_vs_k3"]
        for key, lbl in (("curve_two_sided_alpha0.05", "two_sided"),
                         ("curve_one_sided_alpha0.05", "one_sided")):
            th = d73["theta_for_80pct_power__" + key]
            mde["MDE_" + lbl + "_theta_for_80pct_power"] = th
            mde["MDE_" + lbl + "_fold_change_for_80pct_power_lognormal"] = (
                theta_to_fold(th, sl) if th is not None else None)
        out["MDE_headline"] = mde

    # ------------------------------------------------------- multiplicity
    def holm(ps):
        idx = sorted(range(len(ps)), key=lambda i: ps[i])
        out = [0.0] * len(ps)
        m = len(ps)
        run = 0.0
        for r, i in enumerate(idx):
            run = max(run, min(1.0, (m - r) * ps[i]))
            out[i] = run
        return out

    mult = {"per_panel_over_the_4_labellings": OrderedDict(), "global": {}}
    allp = []
    for pn, p in out["panels"].items():
        keys = [k for k, s in p["tests"].items()
                if "result" not in s and "piN_piS" in s]
        ps = [p["tests"][k]["piN_piS"]["mannwhitney_exact"]["p_two_sided_exact"]
              for k in keys]
        h = holm(ps)
        mult["per_panel_over_the_4_labellings"][pn] = [
            {"labelling": k, "p_raw_two_sided_exact": a,
             "p_holm_within_panel": b, "p_bonferroni_within_panel": min(1.0, len(ps) * a)}
            for k, a, b in zip(keys, ps, h)]
        allp += [(pn, k, a) for k, a in zip(keys, ps)]
    hg = holm([a for _, _, a in allp])
    mult["global"] = {
        "n_piN_piS_group_tests_run": len(allp),
        "tests": [{"panel": pn, "labelling": k, "p_raw": a, "p_holm_global": b}
                  for (pn, k, a), b in zip(allp, hg)],
        "n_raw_p_below_0.05": sum(1 for _, _, a in allp if a < 0.05),
        "n_holm_global_below_0.05": sum(1 for b in hg if b < 0.05),
        "note": ("The pre-registered test is the `disease` labelling on panels P1 "
                 "and P2 only.  The other labellings and panels are the "
                 "sensitivity runs the brief requires; they are listed here so "
                 "that any reader can apply the multiplicity correction "
                 "themselves rather than being shown a single small p out of "
                 "context."),
    }
    out["multiplicity"] = mult

    with open(os.path.join(outdir, "stats_nonparametric.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print("wrote", os.path.join(outdir, "stats_nonparametric.json"))

    write_tables(out, outdir)
    cross_check(out, outdir)


def stratified_perm(rws, labels, strat_key):
    """
    Exact permutation test of the pathogenic/environmental labels permuted
    ONLY WITHIN strata (host class, or rodent/non-rodent).  Statistic: the
    difference of group medians of piN/piS.  Complete enumeration of the
    stratified label space.
    """
    strata = OrderedDict()
    for r in rws:
        strata.setdefault(strat_key(r), []).append(r)
    info = OrderedDict()
    total_space = 1
    for s, members in strata.items():
        npath = sum(1 for r in members if r["category"] == "pathogenic")
        info[s] = {"n": len(members), "n_pathogenic": npath,
                   "n_environmental": len(members) - npath,
                   "species": [r["abbrev"] for r in members],
                   "n_label_arrangements": math.comb(len(members), npath)}
        total_space *= info[s]["n_label_arrangements"]
    vals = {r["abbrev"]: fnum(r["piN_piS"]) for r in rws}
    obs_p = [vals[r["abbrev"]] for r in rws if r["category"] == "pathogenic"]
    obs_e = [vals[r["abbrev"]] for r in rws if r["category"] == "environmental"]
    res = {"strata": info, "size_of_stratified_permutation_space": total_space}
    if not obs_p or not obs_e:
        res["result"] = "UNDEFINED: one label group is empty"
        return res
    obs = median(obs_e) - median(obs_p)
    # enumerate: for each stratum choose which members get the pathogenic label
    choice_sets = []
    for s, members in strata.items():
        npath = info[s]["n_pathogenic"]
        choice_sets.append([set(c) for c in
                            itertools.combinations(range(len(members)), npath)])
    nulls = []
    for combo in itertools.product(*choice_sets):
        pv, ev = [], []
        for (s, members), chosen in zip(strata.items(), combo):
            for i, r in enumerate(members):
                (pv if i in chosen else ev).append(vals[r["abbrev"]])
        if not pv or not ev:
            continue
        nulls.append(median(ev) - median(pv))
    ge = sum(1 for v in nulls if abs(v) >= abs(obs) - 1e-15)
    res.update({
        "statistic": "median(environmental piN/piS) - median(pathogenic piN/piS)",
        "observed": obs,
        "n_arrangements_enumerated": len(nulls),
        "p_two_sided_exact": ge / len(nulls),
        "p_one_sided_exact_env_gt_path":
            sum(1 for v in nulls if v >= obs - 1e-15) / len(nulls),
        "min_attainable_p_two_sided": 1.0 / len(nulls) if len(nulls) else None,
        "note": ("The permutation space collapses to %d arrangements because the "
                 "labels are almost perfectly nested inside the strata; the "
                 "smallest p this design can return is %s."
                 % (len(nulls),
                    ("%.4f" % (1.0 / len(nulls))) if len(nulls) else "n/a")),
    })
    return res


def _parse_newick(newick):
    i = [0]

    def parse():
        if newick[i[0]] == "(":
            i[0] += 1
            kids = []
            while True:
                kids.append(parse())
                if newick[i[0]] == ",":
                    i[0] += 1
                    continue
                if newick[i[0]] == ")":
                    i[0] += 1
                    break
            while i[0] < len(newick) and newick[i[0]] not in ",);":
                i[0] += 1
            return ("int", kids)
        j = i[0]
        while newick[j] not in ":,);":
            j += 1
        name = newick[i[0]:j]
        i[0] = j
        while i[0] < len(newick) and newick[i[0]] not in ",);":
            i[0] += 1
        return ("tip", name)

    return parse()


def sister_pairs(rws, labels, newick):
    """
    Every node of the ML tree whose two descendant clades contain, between the
    species that HAVE a usable piN/piS value in this panel, exactly one
    pathogenic species on one side and exactly one environmental species on
    the other.  These are the phylogenetically paired pathogenic/environmental
    contrasts; a sign test on them is the nonparametric phylogenetic control.
    """
    root = _parse_newick(newick)
    have = {r["abbrev"]: r for r in rws}

    nodes = []

    def walk(node):
        if node[0] == "tip":
            s = {node[1]}
            nodes_local = []
            return s, nodes_local
        childsets = []
        for k in node[1]:
            s, _ = walk(k)
            childsets.append(s)
        allset = set().union(*childsets)
        if len(childsets) == 2:
            a, b = childsets
            ap = [t for t in a if t in have and have[t]["category"] == "pathogenic"]
            ae = [t for t in a if t in have and have[t]["category"] == "environmental"]
            bp = [t for t in b if t in have and have[t]["category"] == "pathogenic"]
            be = [t for t in b if t in have and have[t]["category"] == "environmental"]
            if (len(ap) == 1 and len(ae) == 0 and len(bp) == 0 and len(be) == 1) or \
               (len(ap) == 0 and len(ae) == 1 and len(bp) == 1 and len(be) == 0):
                p = (ap + bp)[0]
                e = (ae + be)[0]
                nodes.append((p, e, sorted(a), sorted(b)))
        return allset, nodes

    walk(root)
    out = {"n_pairs_found": len(nodes), "pairs": []}
    signs = []
    for p, e, a, b in nodes:
        vp = fnum(have[p]["piN_piS"])
        ve = fnum(have[e]["piN_piS"])
        out["pairs"].append({
            "pathogenic": p, "pathogenic_piN_piS": vp,
            "environmental": e, "environmental_piN_piS": ve,
            "clade_A_tips": a, "clade_B_tips": b,
            "sign": "pathogenic LOWER" if vp < ve else
                    "pathogenic HIGHER" if vp > ve else "tie",
        })
        if vp != ve:
            signs.append(1 if vp < ve else -1)
    k = len(signs)
    if k == 0:
        out["sign_test"] = {"result": "no usable pair; no test possible"}
    else:
        pos = sum(1 for s in signs if s == 1)
        # exact two-sided sign test
        tot = 2 ** k
        cnt = 0
        for j in range(k + 1):
            if min(j, k - j) <= min(pos, k - pos):
                cnt += math.comb(k, j)
        out["sign_test"] = {
            "n_pairs_used": k,
            "n_with_pathogenic_lower": pos,
            "p_two_sided_exact": cnt / tot,
            "min_attainable_p_two_sided": 2.0 / tot if k >= 1 else None,
            "note": ("With %d independent sister contrast(s) the smallest "
                     "attainable two-sided exact sign-test p is %.4f, so this "
                     "control is %s."
                     % (k, 2.0 / tot,
                        "uninformative by construction" if 2.0 / tot > 0.05
                        else "capable of rejecting")),
        }
    return out


_PERM10 = {}


def spearman_exact(x, y):
    """Spearman rho with tie-corrected midranks, and an exact two-sided
    permutation p-value by complete enumeration of the n! relabellings when
    n <= 9, or of 10! = 3628800 when n == 10 (vectorised)."""
    n = len(x)
    rx = midranks(x)
    ry = midranks(y)

    def rho(a, b):
        ma = sum(a) / n
        mb = sum(b) / n
        sa = math.sqrt(sum((v - ma) ** 2 for v in a))
        sb = math.sqrt(sum((v - mb) ** 2 for v in b))
        if sa == 0 or sb == 0:
            return None
        return sum((a[i] - ma) * (b[i] - mb) for i in range(n)) / (sa * sb)

    obs = rho(rx, ry)
    if obs is None:
        return {"n": n, "rho": None, "result": "a variable is constant"}
    import numpy as np
    if n > 10:
        # n! too large to enumerate: Monte-Carlo permutation with a fixed seed
        NMC = 1000000
        rg = np.random.default_rng([MASTER_SEED, 77, n])
        a = np.asarray(rx, dtype=float)
        b = np.asarray(ry, dtype=float)
        a = (a - a.mean()) / np.sqrt(((a - a.mean()) ** 2).sum())
        b = (b - b.mean()) / np.sqrt(((b - b.mean()) ** 2).sum())
        P = np.argsort(rg.random((NMC, n)), axis=1)
        rhos = (a[P] * b[None, :]).sum(axis=1)
        ge = int((np.abs(rhos) >= abs(obs) - 1e-12).sum())
        return {"n": n, "rho": obs,
                "n_permutations_sampled": NMC,
                "p_two_sided_monte_carlo": (ge + 1) / (NMC + 1),
                "p_note": ("n! is too large to enumerate; %d random "
                           "permutations, rng seeded [%d, 77, %d]; the "
                           "reported p uses the (ge+1)/(N+1) convention."
                           % (NMC, MASTER_SEED, n)),
                "has_ties_x": len(set(x)) != n, "has_ties_y": len(set(y)) != n}
    key = n
    if key not in _PERM10:
        _PERM10[key] = np.array(list(itertools.permutations(range(n))),
                                dtype=np.int8)
    P = _PERM10[key]
    a = np.asarray(rx, dtype=float)
    b = np.asarray(ry, dtype=float)
    a = (a - a.mean()) / np.sqrt(((a - a.mean()) ** 2).sum())
    b = (b - b.mean()) / np.sqrt(((b - b.mean()) ** 2).sum())
    rhos = (a[P] * b[None, :]).sum(axis=1)
    ge = int((np.abs(rhos) >= abs(obs) - 1e-12).sum())
    tot = int(P.shape[0])
    return {"n": n, "rho": obs,
            "n_permutations_enumerated": tot,
            "p_two_sided_exact": ge / tot,
            "has_ties_x": len(set(x)) != n, "has_ties_y": len(set(y)) != n}


def spearman_block(rws, panels, pname):
    """
    Spearman correlations of per-species piN/piS against three candidate
    confounders.  Sampling depth is always taken from the species' FULL
    single-population set (panel P1), because on an equal-n panel every
    n_sequences is identical by construction.
    """
    # depth/metadata come from the FULL-level panel of the SAME sequence set
    seqset = panels[pname]["seqset"]
    full = {}
    for q in panels.values():
        if q["seqset"] == seqset and q["n_level"] == "full":
            full = {r["abbrev"]: r for r in q["rows"]}
            break
    res = OrderedDict()
    abb = [r["abbrev"] for r in rws]
    y = [fnum(r["piN_piS"]) for r in rws]
    covars = OrderedDict()
    depth, hum, ctry, pis = [], [], [], []
    ok_depth = True
    for r in rws:
        fr = full.get(r["abbrev"])
        if fr is None:
            ok_depth = False
            break
        ns = fnum(fr["n_sequences"])
        nh = fnum(fr["n_human"])
        nc = fnum(fr["n_countries"])
        if None in (ns, nh, nc):
            ok_depth = False
            break
        depth.append(ns)
        hum.append(nh / ns)
        ctry.append(nc)
        pis.append(fnum(r["piS"]))
    if ok_depth:
        covars["n_sequences_in_full_single_population_set"] = depth
        covars["fraction_human_derived_in_full_set"] = hum
        covars["n_countries_in_full_set"] = ctry
        covars["piS_same_panel"] = pis
    for cname, cv in covars.items():
        res[cname] = spearman_exact(cv, y)
        res[cname]["species_order"] = abb
        res[cname]["covariate_values"] = cv
    res["__note"] = (
        "A significant positive correlation with sampling depth or with the "
        "human-derived fraction would mean the per-species piN/piS is driven by "
        "study effort or by clinical/outbreak clonality rather than by selection. "
        "Sampling depth is taken from the full single-population set even on "
        "equal-n panels.")
    return res


def nw_hps_clade(newick):
    """Return the tip set of the smallest clade containing SNV in the given
    newick string.  Pure-python parse; no dendropy/ete3 needed."""
    i = 0

    def parse():
        nonlocal i
        if newick[i] == "(":
            i += 1
            kids = []
            while True:
                kids.append(parse())
                if newick[i] == ",":
                    i += 1
                    continue
                if newick[i] == ")":
                    i += 1
                    break
            # skip label/branch length
            while i < len(newick) and newick[i] not in ",)" and newick[i] != ";":
                i += 1
            return ("int", kids)
        j = i
        while newick[j] not in ":,);":
            j += 1
        name = newick[i:j]
        i = j
        while i < len(newick) and newick[i] not in ",);":
            i += 1
        return ("tip", name)

    root = parse()

    best = None

    def walk(node):
        nonlocal best
        if node[0] == "tip":
            return {node[1]}
        s = set()
        for k in node[1]:
            s |= walk(k)
        if "SNV" in s and (best is None or len(s) < len(best)) and len(s) > 1:
            best = set(s)
        return s

    walk(root)
    # the smallest clade containing SNV is (SNV,(ANDV,LANV),CHOV)-style; walk up
    # one step is not needed: we want the New World HPS clade = smallest clade
    # containing SNV, ANDV and CHOV.
    best2 = None

    def walk2(node):
        nonlocal best2
        if node[0] == "tip":
            return {node[1]}
        s = set()
        for k in node[1]:
            s |= walk2(k)
        if {"SNV", "ANDV", "CHOV"} <= s and (best2 is None or len(s) < len(best2)):
            best2 = set(s)
        return s

    walk2(root)
    return best2


def write_tables(out, outdir):
    """Flat TSV of every test actually run, for the manuscript table."""
    path = os.path.join(outdir, "tests_summary.tsv")
    cols = ["panel", "seqset", "n_level", "labelling", "variable",
            "n_group1", "n_group2", "group1_species", "group2_species",
            "median_group1", "median_group2", "mean_group1", "mean_group2",
            "direction_median", "direction_mean",
            "W_ranksum_group1", "U_group1_over_group2", "theta_P(g2>g1)",
            "p_two_sided_exact", "p_one_sided_exact_g1_lt_g2",
            "min_attainable_p_two_sided", "min_attainable_p_one_sided",
            "n_permutations", "HL_shift_g2_minus_g1", "HL_ci_lo", "HL_ci_hi",
            "HL_ci_level", "perm_diff_median_p_two", "perm_diff_mean_p_two",
            "perm_ratio_mean_p_two", "significant_at_0.05"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for pname, p in out["panels"].items():
            for lkey, sub in p["tests"].items():
                if "result" in sub:
                    continue
                lname = lkey.split(" (")[0]
                for var, t in sub.items():
                    if not isinstance(t, dict) or "mannwhitney_exact" not in t:
                        continue
                    mw = t["mannwhitney_exact"]
                    hl = t["hodges_lehmann"]
                    pm = t["label_permutation_exact"]
                    row = [pname, p["definition"]["seqset"], p["definition"]["n_level"],
                           lname, var,
                           mw["n_x"], mw["n_y"],
                           ",".join(sub["group1_species"]), ",".join(sub["group2_species"]),
                           t["group1"]["median"], t["group2"]["median"],
                           t["group1"]["mean"], t["group2"]["mean"],
                           t["direction_median"], t["direction_mean"],
                           mw["rank_sum_x_W"], mw["U_x_over_y"], mw["theta_hat_P(y>x)"],
                           mw["p_two_sided_exact"], mw["p_one_sided_exact_x_lt_y"],
                           mw["min_attainable_p_two_sided"], mw["min_attainable_p_one_sided"],
                           mw["n_permutations_enumerated"],
                           hl["HL_shift_y_minus_x"], hl["ci_lo"], hl["ci_hi"],
                           hl["ci_level_achieved"],
                           pm["diff_median"]["p_two_sided_exact"],
                           pm["diff_mean"]["p_two_sided_exact"],
                           pm["ratio_mean"]["p_two_sided_exact"],
                           "YES" if mw["p_two_sided_exact"] < 0.05 else "no"]
                    fh.write("\t".join("" if v is None else
                                       ("%.6g" % v if isinstance(v, float) else str(v))
                                       for v in row) + "\n")
    print("wrote", path)


def cross_check(out, outdir):
    """Independent cross-check of the headline exact p-values against scipy,
    if scipy is importable.  Written to its own file; never used to produce a
    reported number."""
    res = {"scipy_available": False}
    try:
        import scipy
        from scipy import stats as sps
        res["scipy_available"] = True
        res["scipy_version"] = scipy.__version__
        checks = []
        for pname, p in out["panels"].items():
            for lkey, sub in p["tests"].items():
                if "result" in sub:
                    continue
                t = sub.get("piN_piS")
                if not t:
                    continue
                g1 = [s["piN_piS"] for s in p["species"]
                      if s["abbrev"] in sub["group1_species"]]
                g2 = [s["piN_piS"] for s in p["species"]
                      if s["abbrev"] in sub["group2_species"]]
                try:
                    r = sps.mannwhitneyu(g1, g2, alternative="two-sided",
                                         method="exact")
                    sp_p = float(r.pvalue)
                    sp_u = float(r.statistic)
                except Exception as e:  # noqa
                    sp_p, sp_u = None, str(e)
                mine = t["mannwhitney_exact"]
                checks.append({
                    "panel": pname, "labelling": lkey,
                    "my_U": mine["U_x_over_y"], "scipy_U": sp_u,
                    "my_p_two": mine["p_two_sided_exact"], "scipy_p_two": sp_p,
                    "agree": (sp_p is not None
                              and abs(sp_p - mine["p_two_sided_exact"]) < 1e-12
                              and abs(sp_u - mine["U_x_over_y"]) < 1e-9),
                })
        res["checks"] = checks
        res["n_checks"] = len(checks)
        res["n_agree"] = sum(1 for c in checks if c["agree"])
    except ImportError as e:
        res["error"] = str(e)
    with open(os.path.join(outdir, "crosscheck_scipy.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print("wrote", os.path.join(outdir, "crosscheck_scipy.json"),
          "agree %s/%s" % (res.get("n_agree"), res.get("n_checks")))


if __name__ == "__main__":
    main()
