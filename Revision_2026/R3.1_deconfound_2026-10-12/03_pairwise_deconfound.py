#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""R3.1 + R3.3c on the pairwise table (the manuscript's own level of analysis).

Four things, all on the SAME 231 pairs from
R2.2/results/pairwise_three_estimators.tsv (manuscript pn/ps, NG86, YN00):

 A. HEAD-TO-HEAD PARTITIONS. Relabel the identical 231 pairs by HOST
    (rodent/non-rodent) instead of by DISEASE (pathogenic/environmental) and
    compare Kruskal-Wallis H and group means. Same data, same n, so H is
    directly comparable. If host separates the groups as well as or better than
    disease, the manuscript's attribution is not supported.

 B. RESTRICTED (within-host-stratum) EXACT PERMUTATION TEST. The manuscript
    permutes disease labels over all 22 species, which moves labels across host
    classes and so cannot separate the two. Permuting disease labels ONLY
    within the rodent-borne stratum holds host structure fixed; since all 10
    pathogenic species are rodent-borne and there are 12 rodent-borne species,
    the restricted permutation set is exactly C(12,10) = 66 assignments, so the
    test can be enumerated exhaustively. Minimum attainable p = 1/66 = 0.0152.

 C. DIVERGENCE-CONTROLLED regression: omega ~ dS + group, because the
    environmental-environmental pairs are also the deepest (R3.2's saturation
    point). Significance of the group term by the restricted permutation of B.

 D. R3.3c SENSITIVITY: TULV reassigned / dropped, plus a leave-one-species-out
    jackknife over all 22 species, asking in each case whether the ordering
    path-path < cross < env-env survives.
"""
import itertools
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

A = ("/PATH/TO/"
     "Project_Hantavirus/Reviews_ans/analysis/R3.1_deconfound_2026-10-12")

PATHOGENIC = ["HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV",
              "CHOV", "LANV"]
RODENT = PATHOGENIC + ["PHV", "TULV"]
ALL22 = RODENT + ["TPMV", "MJNV", "NVAV", "ASAV", "OXBV", "RKPV", "SWSV",
                  "CBNV", "QZNV", "JJUV"]

df = pd.read_csv(f"{A}/inputs/pairwise_three_estimators.tsv", sep="\t")
assert len(df) == 231, len(df)
EST = {"ms_omega": "manuscript pn/ps (uncorrected)",
       "ng_omega": "NG86 (JC-corrected)",
       "yn_omega": "YN00 (ts/tv + codon usage)"}

out = []
def say(s=""):
    out.append(s)
    print(s)


def pair_group(a, b, members):
    """3-level pair label given the set of 'focal' species."""
    ia, ib = a in members, b in members
    return "in-in" if (ia and ib) else ("out-out" if not (ia or ib) else "cross")


def summarise(col, groups, keep=None):
    """mean/median/n per group + Kruskal-Wallis over the 3 groups."""
    g = pd.Series(groups, index=df.index)
    d = pd.DataFrame({"v": df[col], "g": g})
    if keep is not None:
        d = d[keep]
    d = d.dropna(subset=["v"])
    res = {}
    for name in ("in-in", "cross", "out-out"):
        v = d.loc[d.g == name, "v"].values
        res[name] = dict(n=int(len(v)),
                         mean=float(np.mean(v)) if len(v) else float("nan"),
                         median=float(np.median(v)) if len(v) else float("nan"))
    arrs = [d.loc[d.g == n, "v"].values for n in ("in-in", "cross", "out-out")]
    if all(len(a) > 0 for a in arrs):
        H, p = stats.kruskal(*arrs)
    else:
        H, p = float("nan"), float("nan")
    res["H"] = float(H)
    res["p"] = float(p)
    res["n_total"] = int(len(d))
    res["ordering_preserved"] = bool(
        res["in-in"]["mean"] < res["cross"]["mean"] < res["out-out"]["mean"])
    return res


say("# R3.1 / R3.3c — pairwise-level deconfounding and sensitivity")
say()
say("Input: the same 231 pairs and the same three estimators used for R3.2")
say("(`R2.2_phylo_dnds_2026-10-12/results/pairwise_three_estimators.tsv`).")
say("`in-in` = both species in the focal set, `out-out` = neither.")
say()

# ---------------------------------------------------------------- A -----------
say("## A. Head-to-head: DISEASE partition vs HOST partition, identical pairs")
say()
dis = [pair_group(r.species_A, r.species_B, set(PATHOGENIC))
       for r in df.itertuples()]
hos = [pair_group(r.species_A, r.species_B, set(RODENT))
       for r in df.itertuples()]
say("Partition sizes (pairs): disease in-in/cross/out-out = "
    f"{dis.count('in-in')}/{dis.count('cross')}/{dis.count('out-out')}; "
    f"host = {hos.count('in-in')}/{hos.count('cross')}/{hos.count('out-out')}")
say()
headA = {}
for col, lab in EST.items():
    rd = summarise(col, dis)
    rh = summarise(col, hos)
    headA[col] = dict(disease=rd, host=rh)
    say(f"### {lab} — `{col}`")
    say()
    say("| partition | in-in mean | cross mean | out-out mean | n | KW H | p | ordering in<cross<out |")
    say("|---|---|---|---|---|---|---|---|")
    for nm, r in (("disease (pathogenic)", rd), ("host (rodent-borne)", rh)):
        say(f"| {nm} | {r['in-in']['mean']:.4f} | {r['cross']['mean']:.4f} | "
            f"{r['out-out']['mean']:.4f} | {r['n_total']} | {r['H']:.3f} | "
            f"{r['p']:.3e} | {r['ordering_preserved']} |")
    better = "HOST" if rh["H"] > rd["H"] else "DISEASE"
    say()
    say(f"Larger Kruskal-Wallis H (same pairs, same n = {rd['n_total']}): "
        f"**{better}** (disease H = {rd['H']:.3f}, host H = {rh['H']:.3f})")
    say()

# ---------------------------------------------------------------- B -----------
say("## B. Restricted exact permutation: disease labels permuted WITHIN the "
    "rodent-borne stratum only")
say()
say("All 10 pathogenic species are rodent-borne and 12 species are "
    "rodent-borne, so holding the host structure fixed leaves exactly "
    "C(12,10) = 66 possible disease labellings. The test below enumerates all "
    "66; the observed labelling is one of them, so the smallest attainable "
    "one-sided p is 1/66 = 0.0152.")
say()
assignments = list(itertools.combinations(RODENT, 10))
assert len(assignments) == 66
obs_set = tuple(PATHOGENIC)
headB = {}
for col, lab in EST.items():
    stat_obs = None
    stats_all = []
    for asg in assignments:
        g = [pair_group(r.species_A, r.species_B, set(asg))
             for r in df.itertuples()]
        r = summarise(col, g)
        # statistic: how much LOWER in-in is than out-out (manuscript's claim)
        s = r["out-out"]["mean"] - r["in-in"]["mean"]
        stats_all.append(s)
        if set(asg) == set(obs_set):
            stat_obs = s
    stats_all = np.array(stats_all)
    assert stat_obs is not None
    # one-sided: how many restricted labellings give a separation at least as
    # extreme as the observed one
    ge = int(np.sum(stats_all >= stat_obs - 1e-12))
    p_exact = ge / len(stats_all)
    rank = int(np.sum(stats_all > stat_obs + 1e-12)) + 1
    headB[col] = dict(stat_obs=float(stat_obs), p_exact=p_exact, rank=rank,
                      n_perm=len(stats_all),
                      max=float(stats_all.max()), mean=float(stats_all.mean()))
    say(f"- **{lab}**: observed separation (out-out mean - in-in mean) = "
        f"{stat_obs:.4f}; rank {rank} of 66 restricted labellings "
        f"(1 = most extreme); exact one-sided p = **{p_exact:.4f}** "
        f"(null mean {stats_all.mean():.4f}, max {stats_all.max():.4f})")
say()
say("Interpretation: this p is the probability of a separation this large when "
    "the disease labels are shuffled among rodent-borne species only. Unlike "
    "the manuscript's unrestricted permutation it cannot be driven by host "
    "class, because host class never changes. It is, however, a test with only "
    "66 possible outcomes and the 66 labellings overlap heavily with the "
    "observed one (8 of 10 pathogenic species are shared by any two of them), "
    "so it has little power. A null result here is weak evidence either way, "
    "and we report it as such.")
say()

# ---------------------------------------------------------------- C -----------
say("## C. Divergence-controlled regression (R3.2's saturation point as a "
    "confound in its own right)")
say()
say("The environmental-environmental pairs are also the deepest pairs, so "
    "group and divergence are themselves confounded. OLS of omega on dS plus "
    "group indicators, with the group term tested by the restricted "
    "permutation of section B.")
say()
headC = {}
for col, dscol in (("yn_omega", "yn_dS"), ("ms_omega", "ms_ps")):
    sub = df[[ "species_A", "species_B", col, dscol]].dropna()
    y = sub[col].values
    dS = sub[dscol].values

    def design(members):
        g = [pair_group(a, b, set(members))
             for a, b in zip(sub.species_A, sub.species_B)]
        inin = np.array([x == "in-in" for x in g], float)
        outout = np.array([x == "out-out" for x in g], float)
        return np.column_stack([np.ones_like(dS), dS, inin, outout])

    def fit(members):
        X = design(members)
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        rss = float(resid @ resid)
        return beta, rss

    b_obs, rss_obs = fit(PATHOGENIC)
    # statistic: contrast between in-in and out-out coefficients after dS
    contrast_obs = b_obs[3] - b_obs[2]
    contrasts = []
    for asg in assignments:
        b, _ = fit(list(asg))
        contrasts.append(b[3] - b[2])
    contrasts = np.array(contrasts)
    p_exact = float(np.sum(contrasts >= contrast_obs - 1e-12) / len(contrasts))
    # how much variance does dS alone explain
    X0 = np.column_stack([np.ones_like(dS), dS])
    b0, *_ = np.linalg.lstsq(X0, y, rcond=None)
    rss0 = float(((y - X0 @ b0) ** 2).sum())
    tss = float(((y - y.mean()) ** 2).sum())
    headC[col] = dict(n=int(len(y)), beta_dS=float(b_obs[1]),
                      contrast=float(contrast_obs), p_exact=p_exact,
                      r2_dS_only=1 - rss0 / tss, r2_full=1 - rss_obs / tss)
    say(f"- **{EST[col]}** (n = {len(y)} pairs): slope on {dscol} = "
        f"{b_obs[1]:.4f}; R2 with divergence alone = {1 - rss0 / tss:.3f}, "
        f"with divergence + disease group = {1 - rss_obs / tss:.3f}; "
        f"out-out minus in-in group effect after divergence = "
        f"{contrast_obs:.4f}, restricted-permutation one-sided p = "
        f"**{p_exact:.4f}**")
say()

# ---------------------------------------------------------------- D -----------
say("## D. R3.3c — reassignment sensitivity")
say()
say("The manuscript's own metadata is internally inconsistent for Tula: "
    "`Work1/Genus_Diversity/genus_diversity_manifest.tsv` records TULV's note "
    "as \"Rare/mild human case reports only\" while its `category` column "
    "reads `environmental`. The contradiction in the prose is "
    "therefore present in the data files too.")
say()
scen = {
    "S0 baseline": (PATHOGENIC, None),
    "S1 TULV -> pathogenic": (PATHOGENIC + ["TULV"], None),
    "S2 TULV dropped": (PATHOGENIC, "TULV"),
    "S3 TULV + PHV -> pathogenic (= host partition)": (RODENT, None),
}
headD = {}
for col, lab in EST.items():
    say(f"### {lab} — `{col}`")
    say()
    say("| scenario | in-in mean | cross mean | out-out mean | n | KW H | p | ordering preserved |")
    say("|---|---|---|---|---|---|---|---|")
    for nm, (members, drop) in scen.items():
        keep = None
        if drop:
            keep = ~((df.species_A == drop) | (df.species_B == drop))
        g = [pair_group(r.species_A, r.species_B, set(members))
             for r in df.itertuples()]
        r = summarise(col, g, keep=keep)
        headD[(col, nm)] = r
        say(f"| {nm} | {r['in-in']['mean']:.4f} | {r['cross']['mean']:.4f} | "
            f"{r['out-out']['mean']:.4f} | {r['n_total']} | {r['H']:.3f} | "
            f"{r['p']:.3e} | {r['ordering_preserved']} |")
    say()

say("### Leave-one-species-out jackknife (all 22 deletions, baseline labels)")
say()
jack = {}
for col, lab in EST.items():
    rows = []
    for sp in ALL22:
        keep = ~((df.species_A == sp) | (df.species_B == sp))
        g = [pair_group(r.species_A, r.species_B, set(PATHOGENIC))
             for r in df.itertuples()]
        r = summarise(col, g, keep=keep)
        rows.append((sp, r))
    bad = [sp for sp, r in rows if not r["ordering_preserved"]]
    ps = [r["p"] for _, r in rows]
    jack[col] = dict(broken=bad, p_max=float(np.nanmax(ps)),
                     p_min=float(np.nanmin(ps)))
    say(f"- **{lab}**: ordering path-path < cross < env-env broken by "
        f"{len(bad)} of 22 single-species deletions"
        + (f" ({', '.join(bad)})" if bad else "")
        + f"; Kruskal-Wallis p across the 22 jackknife replicates ranges "
          f"{np.nanmin(ps):.2e} to {np.nanmax(ps):.2e}")
say()

with open(f"{A}/results/03_pairwise_deconfound.md", "w") as fh:
    fh.write("\n".join(out) + "\n")
with open(f"{A}/results/03_pairwise_deconfound.json", "w") as fh:
    json.dump({"A": headA, "B": headB, "C": headC,
               "D": {f"{k[0]}|{k[1]}": v for k, v in headD.items()},
               "jackknife": jack}, fh, indent=1)
print("\nwrote results/03_pairwise_deconfound.{md,json}")
