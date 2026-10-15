"""Diagnostics for the phylogenetic treatment: the host3 singleton-group artefact,
species jackknife of the disease contrast, and the AICc model ranking.

Usage: python3 -I 24_diagnostics.py <ANALYSIS_ROOT> <LABELS_TSV> <TREE_NWK> <OUTDIR>
"""
import csv
import importlib.util
import json
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


def build(root, labels_tsv, tree_nwk, seqset, level):
    labels = ps.load_labels(labels_tsv)
    rows = ps.pick(ps.load_table(os.path.join(root, "results",
                                              "pinps_per_species.tsv")),
                   seqset, level, True)
    species = sorted(r["abbrev"] for r in rows)
    byab = {r["abbrev"]: r for r in rows}
    y = np.array([float(byab[s]["piN_piS"]) for s in species])
    tree = pt.parse_newick(open(tree_nwk).read())
    return labels, species, y, tree


def vcv_for(tree, species):
    pr = pt.prune(tree, set(species))
    pr.parent = None
    adj, names = pt.adjacency(pr)
    return pt.vcv_rooted_at_node(adj, names, species, 0), pr


def main():
    root, labels_tsv, tree_nwk, outdir = sys.argv[1:5]
    labels, species, y, tree = build(root, labels_tsv, tree_nwk, "clust95d", "10")
    C, pruned = vcv_for(tree, species)
    out = dict(species=species, y=dict(zip(species, [float(v) for v in y])))

    # ---------- host3 singleton-group artefact
    X, cn, g = ps.design(species, "host3", labels)
    f = ps.gls(y, X, C)
    fitted = X @ f["beta"]
    out["host3"] = dict(
        colnames=cn, beta=[float(b) for b in f["beta"]],
        group_counts={h: int(sum(1 for s in species if labels[s]["host_class"] == h))
                      for h in sorted({labels[s]["host_class"] for s in species})},
        fitted=dict(zip(species, [float(v) for v in fitted])),
        residual=dict(zip(species, [float(v) for v in (y - fitted)])),
        note=("each non-rodent host class contains exactly one species, so the two "
              "extra parameters fit those species' values exactly; the fit gains "
              "nothing about host biology"),
    )

    # ---------- AICc ranking, piN/piS, BM
    rank = []
    for lab in ["intercept", "disease", "disease_tulv_moved", "host2", "host3",
                "disease+host2", "cladeNW"]:
        X, cn, g = ps.design(species, lab, labels)
        fl = ps.gls(y, X, C)
        k = X.shape[1] + 1
        aic, aicc = ps.aic_set(fl["logL"], k, len(y))
        rank.append(dict(model=lab, k_params=k, logL=fl["logL"], AIC=aic, AICc=aicc,
                         colnames=cn))
    a0 = [r for r in rank if r["model"] == "intercept"][0]["AICc"]
    for r in rank:
        r["dAICc_vs_intercept"] = r["AICc"] - a0
    rank.sort(key=lambda r: r["AICc"])
    out["aicc_ranking_BM_piN_piS"] = rank

    # ---------- GLS weights: which species actually carry each contrast
    out["contrast_weights"] = {}
    Vi = np.linalg.inv(C)
    for lab in ["disease", "host2", "cladeNW"]:
        X, cn, gg = ps.design(species, lab, labels)
        P = np.linalg.inv(X.T @ Vi @ X) @ X.T @ Vi
        Po = np.linalg.inv(X.T @ X) @ X.T
        out["contrast_weights"][lab] = dict(
            gls=dict(zip(species, [float(v) for v in P[1]])),
            ols=dict(zip(species, [float(v) for v in Po[1]])),
            beta_gls=float(P[1] @ y), beta_ols=float(Po[1] @ y),
            note=("GLS weights are the second row of (X'V^-1X)^-1 X'V^-1; they sum to "
                  "+1 over the labelled group and -1 over the rest. They say which "
                  "species the phylogenetic contrast is actually made of."))

    # ---------- species jackknife of the disease contrast
    jk = []
    for drop in species:
        keep = [s for s in species if s != drop]
        yk = np.array([out["y"][s] for s in keep])
        Ck, _ = vcv_for(tree, keep)
        X, cn, gk = ps.design(keep, "disease", labels)
        if gk.sum() in (0, len(gk)):
            jk.append(dict(dropped=drop, status="group empty"))
            continue
        fk = ps.gls(yk, X, Ck)
        a, b = yk[gk == 1], yk[gk == 0]
        try:
            mw = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
            U, pU = float(mw.statistic), float(mw.pvalue)
        except Exception:
            U, pU = float("nan"), float("nan")
        jk.append(dict(dropped=drop, k=len(keep), k_env=int(gk.sum()),
                       beta=float(fk["beta"][1]), se=float(fk["se"][1]),
                       t=float(fk["t"][1]), p=float(fk["p"][1]),
                       mean_env=float(a.mean()), mean_path=float(b.mean()),
                       mannwhitney_U=U, mannwhitney_p=pU))
    out["jackknife_disease"] = jk

    # ---------- what is attainable with 2 environmental species (MJNV dropped)
    out["note_two_env"] = ("with 2 environmental vs 7 pathogenic species the "
                           "smallest attainable two-sided exact Mann-Whitney p is "
                           "2/C(9,2) = %.4f" % (2.0 / 36))

    with open(os.path.join(outdir, "diagnostics.json"), "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    with open(os.path.join(outdir, "aicc_ranking.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["model", "k_params", "logL", "AIC", "AICc", "dAICc_vs_intercept"])
        for r in rank:
            w.writerow([r["model"], r["k_params"], "%.4f" % r["logL"],
                        "%.3f" % r["AIC"], "%.3f" % r["AICc"],
                        "%+.3f" % r["dAICc_vs_intercept"]])
    print(json.dumps(out["host3"], indent=1))
    print(json.dumps(rank, indent=1))
    print(json.dumps(jk, indent=1))
    print(out["note_two_env"])


if __name__ == "__main__":
    main()
