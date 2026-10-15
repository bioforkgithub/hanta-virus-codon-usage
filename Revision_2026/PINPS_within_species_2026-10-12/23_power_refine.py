"""Refined power / minimum-detectable-effect calculation for the phylogenetic
treatment of within-species piN/piS, plus confidence intervals on every group
coefficient and a calibration check of the label-randomisation test.

Usage: python3 -I 23_power_refine.py <ANALYSIS_ROOT> <LABELS_TSV> <TREE_NWK> <OUTDIR>
Seed 20261012.
"""
import csv
import importlib.util
import itertools
import json
import math
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("pt", os.path.join(HERE, "20_phylo_tree.py"))
pt = importlib.util.module_from_spec(_s)
_s.loader.exec_module(pt)
_s2 = importlib.util.spec_from_file_location("ps", os.path.join(HERE, "21_phylo_stats.py"))
ps = importlib.util.module_from_spec(_s2)
_s2.loader.exec_module(ps)

SEED = 20261012
NSIM = 20000


def main():
    root, labels_tsv, tree_nwk, outdir = sys.argv[1:5]
    labels = ps.load_labels(labels_tsv)
    rows = ps.pick(ps.load_table(os.path.join(root, "results", "pinps_per_species.tsv")),
                   "clust95d", "10", True)
    species = sorted(r["abbrev"] for r in rows)
    byab = {r["abbrev"]: r for r in rows}
    y = np.array([float(byab[s]["piN_piS"]) for s in species])
    tree = pt.parse_newick(open(tree_nwk).read())
    pruned = pt.prune(tree, set(species))
    pruned.parent = None
    adj, names = pt.adjacency(pruned)
    # BM-PGLS point estimates and standard errors are invariant to root placement
    # (rerooting adds a constant to every element of V, which an intercept absorbs),
    # verified in phylo_stats.json: identical beta/se/t under all three rootings.
    C = pt.vcv_rooted_at_node(adj, names, species, 0)

    out = dict(seed=SEED, species=species, k_species=len(species),
               y_piN_piS={s: float(byab[s]["piN_piS"]) for s in species},
               rooting="as_given (BM estimates are rooting-invariant)")

    # ---- confidence intervals on every group coefficient
    out["coefficients"] = []
    for lab in ["disease", "disease_tulv_moved", "host2", "cladeNW"]:
        X, cn, g = ps.design(species, lab, labels)
        f = ps.gls(y, X, C)
        tc = stats.t.ppf(0.975, f["df"])
        b, se = float(f["beta"][1]), float(f["se"][1])
        a = y[g == 1]
        bb = y[g == 0]
        out["coefficients"].append(dict(
            model=lab, k_group1=int(len(a)), k_group0=int(len(bb)),
            mean_group1=float(a.mean()), mean_group0=float(bb.mean()),
            median_group1=float(np.median(a)), median_group0=float(np.median(bb)),
            beta=b, se=se, t=float(f["t"][1]), p=float(f["p"][1]), df=f["df"],
            ci95=[b - tc * se, b + tc * se], t_crit=float(tc),
            half_width=float(tc * se)))

    # ---- refined power curve
    X, cn, g = ps.design(species, "disease", labels)
    f0 = ps.gls(y, np.ones((len(y), 1)), C)
    sigma2 = f0["sigma2_reml"]
    path_mean = float(y[g == 0].mean())
    grid = [0.0, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225,
            0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.045, 0.05, 0.055,
            0.06, 0.07, 0.08, 0.09, 0.10]
    # add the effect size that matches the BETWEEN-species relative difference
    # (S segment YN00: pathogenic-pathogenic 0.0678 vs environmental-environmental
    #  0.1005, i.e. +48.2% relative) applied to the pathogenic within-species mean
    rel_between = (0.1005 - 0.0678) / 0.0678
    d_between = rel_between * path_mean
    grid.append(round(d_between, 7))
    grid = sorted(set(grid))
    rng = np.random.default_rng([SEED, 3])
    curve = ps.power_curve(C, g, rng, sigma2, grid, nsim=NSIM)
    out["power"] = dict(
        nsim=NSIM, alpha=0.05, sigma2_reml_BM=float(sigma2),
        pathogenic_mean=path_mean,
        environmental_mean=float(y[g == 1].mean()),
        observed_diff=float(y[g == 1].mean() - y[g == 0].mean()),
        between_species_relative_effect=rel_between,
        d_matching_between_species_relative_effect=float(d_between),
        curve=curve)

    # interpolate the 80% power thresholds
    def thresh(key, target=0.8):
        xs = [c["d"] for c in curve]
        ys = [c[key] for c in curve]
        for i in range(1, len(xs)):
            if ys[i - 1] < target <= ys[i]:
                return xs[i - 1] + (target - ys[i - 1]) * (xs[i] - xs[i - 1]) / (
                    ys[i] - ys[i - 1])
        return None

    d80r = thresh("power_randomisation")
    d80m = thresh("power_mannwhitney")
    out["mde"] = dict(
        d80_randomisation=d80r, d80_mannwhitney=d80m,
        ratio80_randomisation=(1 + d80r / path_mean) if d80r else None,
        ratio80_mannwhitney=(1 + d80m / path_mean) if d80m else None,
        power_at_observed=[c for c in curve
                           if abs(c["d"] - out["power"]["observed_diff"]) < 1e-9],
        power_at_between_species_effect=[c for c in curve
                                         if abs(c["d"] - d_between) < 1e-7],
        size_at_d0=[c for c in curve if c["d"] == 0.0],
    )

    with open(os.path.join(outdir, "power_mde.json"), "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    with open(os.path.join(outdir, "power_curve.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["d_env_minus_path", "ratio_env_over_path_mean",
                    "power_label_randomisation_PGLS_t", "power_exact_mannwhitney"])
        for c in curve:
            w.writerow(["%.6f" % c["d"], "%.3f" % (1 + c["d"] / path_mean),
                        "%.4f" % c["power_randomisation"],
                        "%.4f" % c["power_mannwhitney"]])

    print(json.dumps(out["coefficients"], indent=1))
    print(json.dumps(out["mde"], indent=1))
    print("wrote", os.path.join(outdir, "power_mde.json"),
          os.path.join(outdir, "power_curve.tsv"))


if __name__ == "__main__":
    main()
