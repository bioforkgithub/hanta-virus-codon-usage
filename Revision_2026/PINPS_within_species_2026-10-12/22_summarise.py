"""Print the key tables from phylo_stats.json in readable form.

Usage: python3 -I 22_summarise.py <phylo_stats.json> [dataset_tag]
"""
import json
import sys


def f(x, n=4):
    try:
        if x is None:
            return "NA"
        return ("%." + str(n) + "g") % float(x)
    except Exception:
        return str(x)


def main():
    out = json.load(open(sys.argv[1]))
    only = sys.argv[2] if len(sys.argv) > 2 else None
    for ds in out["datasets"]:
        if only and ds.get("tag") != only:
            continue
        if ds.get("skipped"):
            print("## %s  SKIPPED (k=%s)" % (ds["tag"], ds["k_species"]))
            continue
        print("=" * 100)
        print("## DATASET %s   k_species=%d" % (ds["tag"], ds["k_species"]))
        print(ds["notes"])
        print("pruned newick:", ds["pruned_newick"])
        print("rootings:", ds["rootings"], " midpoint tip-to-tip span:",
              f(ds["midpoint_span"]))
        print("crosstab:", ds["crosstab"])
        print("collinearity:", json.dumps(ds["collinearity"]))
        print("log offset:", f(ds["log_offset"], 6))
        print("\n-- per-species input")
        print("%-6s %-14s %-9s %-6s %5s %10s %10s %10s %s" % (
            "sp", "category", "host", "NWcl", "n", "piN", "piS", "piN/piS", "sat"))
        for s in ds["species"]:
            print("%-6s %-14s %-9s %-6s %5d %10s %10s %10s %s" % (
                s, ds["category"][s], ds["host_class"][s],
                "yes" if ds["NW_HPS_clade"][s] else "-",
                ds["n_sequences_used"][s], f(ds["piN"][s], 4), f(ds["piS"][s], 4),
                f(ds["piN_piS"][s], 4), ds["saturation_flags"][s]))

        print("\n-- phylogenetic signal (intercept-only)")
        print("%-20s %-20s %7s %9s %9s %9s %8s %8s %9s %9s %9s" % (
            "rooting", "response", "lamML", "logL_BM", "logL_lam", "logL_star",
            "LRTlam", "pLRT", "K", "AICcBM", "AICcstar"))
        for r in ds["signal"]:
            print("%-20s %-20s %7s %9s %9s %9s %8s %8s %9s %9s %9s" % (
                r["rooting"], r["response"], f(r["lambda_ML"]), f(r["logL_BM"]),
                f(r["logL_lambda"]), f(r["logL_star"]), f(r["LRT_lambda_vs_0"]),
                f(r["p_LRT_lambda_vs_0"]), f(r["blomberg_K"]), f(r["aicc_BM"]),
                f(r["aicc_star"])))

        print("\n-- PGLS fits (group coefficient)")
        print("%-20s %-20s %-20s %-10s %6s %8s %10s %10s %7s %8s %8s %8s %8s" % (
            "rooting", "response", "model", "err", "lam", "logL", "beta", "se",
            "t", "p", "AICc", "LRT", "pLRT"))
        for r in ds["fits"]:
            if r.get("status"):
                print("%-20s %-20s %-20s  %s" % (r["rooting"], r["response"],
                                                 r["model"], r["status"]))
                continue
            gi = 1 if len(r["beta"]) > 1 else 0
            print("%-20s %-20s %-20s %-10s %6s %8s %10s %10s %7s %8s %8s %8s %8s" % (
                r["rooting"], r["response"], r["model"], r["error_model"],
                f(r["lam"], 3), f(r["logL"], 5), f(r["beta"][gi]), f(r["se"][gi]),
                f(r["t"], 3) if not isinstance(r["t"], list) else f(r["t"][gi], 3),
                f(r["p"][gi], 3), f(r["AICc"], 5),
                f(r.get("LRT_vs_intercept"), 3), f(r.get("p_LRT"), 3)))

        print("\n-- Garland phylogenetic ANOVA")
        print("%-20s %-20s %-20s %10s %10s %10s %s" % (
            "rooting", "response", "model", "F", "p_ols", "p_phylo", "groups"))
        for r in ds["anova"]:
            print("%-20s %-20s %-20s %10s %10s %10s %s" % (
                r["rooting"], r["response"], r["model"], f(r["F_obs"]),
                f(r["p_ols"], 3), f(r["p_phylo"], 3), r["group_sizes"]))

        print("\n-- exact label randomisation on the BM-PGLS t statistic")
        print("%-20s %-20s %-20s %5s %5s %8s %8s %10s %8s %6s %8s %8s" % (
            "rooting", "response", "model", "rest", "nass", "pars", "t_obs",
            "beta", "p_t", "p_dif", "MW_U", "MW_p"))
        for r in ds["randomisation"]:
            print("%-20s %-20s %-20s %5s %5d %8s %8s %10s %8s %6s %8s %8s" % (
                r["rooting"], r["response"], r["model"],
                "Y" if r["restricted"] else "N", r["n_assignments"],
                r["parsimony_obs"], f(r["t_obs"], 3), f(r["beta_obs"]),
                f(r["p_t"], 3), f(r["p_diff"], 3), f(r["mannwhitney_U"], 3),
                f(r["mannwhitney_p_exact"], 3)))

        print("\n-- power / minimum detectable effect (response piN/piS, disease labelling)")
        for pw in ds["power"]:
            print("rooting=%s sigma2_reml=%s implied_tip_sd=%s obs_diff=%s nsim=%d" % (
                pw["rooting"], f(pw["sigma2_reml_BM"]), f(pw["implied_tip_sd"]),
                f(pw["observed_diff_env_minus_path"]), pw["nsim"]))
            print("   %10s %12s %12s" % ("d", "pow_rand", "pow_MW"))
            for c in pw["curve"]:
                print("   %10s %12s %12s" % (f(c["d"], 5),
                                             f(c["power_randomisation"], 3),
                                             f(c["power_mannwhitney"], 3)))


if __name__ == "__main__":
    main()
