"""
Exact permutation test for the across-species pn/ps comparison.

WHY THIS EXISTS
---------------
03_genus_stats_figures.R compares the 231 pairwise pn/ps values across three
categories with a Kruskal-Wallis test and post-hoc Wilcoxon tests. Those tests
assume the observations are independent. They are not. Each of the 22 species
contributes to 21 of the 231 pairs, so a single species with unusual pn/ps
moves 21 values at once. The resulting p-values are anti-conservative by
roughly ten orders of magnitude.

The unit that CAN be permuted is the species label, not the pair. This script
enumerates every possible relabelling and computes exact p-values from it.

    C(22, 10) = 646,646 ways to call 10 of the 22 species pathogenic

That is small enough to enumerate completely, so these are exact p-values, not
Monte-Carlo estimates.

TWO OBSERVATIONS MAKE IT CHEAP
------------------------------
1. Kruskal-Wallis ranks the pooled 231 values. The pooled set is the same under
   every relabelling, so the ranks are constant - only group membership moves.
2. Every relabelling keeps 10 pathogenic and 12 environmental species, so the
   group sizes are fixed at 45 / 120 / 66. H therefore reduces to
   sum(R_g^2 / n_g), and the tie correction is a constant that cancels when
   comparing H across permutations.

So the whole null distribution is a few numpy reductions and runs in seconds.

INPUTS
------
    results/pairwise_pnps_matrix.tsv   (from 02_pairwise_dnds.py)
    genus_diversity_manifest.tsv       (species -> category)

OUTPUT
------
    results/permutation_test.tsv
    plus a summary to stdout

Usage:  python 05_permutation_test.py
"""
import csv
import itertools
from math import comb
from pathlib import Path

import numpy as np
from scipy import stats

BASE = Path(__file__).parent.parent
RESULTS = BASE / "results"

CHUNK = 20_000          # relabellings evaluated per numpy pass


def load():
    rows = list(csv.DictReader(
        open(RESULTS / "pairwise_pnps_matrix.tsv"), delimiter="\t"))
    cat = {r["abbrev"]: r["category"] for r in csv.DictReader(
        open(BASE / "genus_diversity_manifest.tsv"), delimiter="\t")}
    species = sorted(cat)
    idx = {s: i for i, s in enumerate(species)}
    A = np.array([idx[r["species_A"]] for r in rows])
    B = np.array([idx[r["species_B"]] for r in rows])
    V = np.array([float(r["pnps"]) for r in rows])
    labels = np.array([cat[s] == "pathogenic" for s in species])
    assert len(V) == comb(len(species), 2), "matrix is not a complete set of pairs"
    return species, A, B, V, labels


def group_stats(L, A, B, V, RANK):
    """L: (C, n_species) bool array of relabellings.

    Returns the Kruskal-Wallis H proxy and the three group means for each row.
    Group code: 2 = pathogenic-pathogenic, 1 = cross, 0 = environmental-environmental.
    """
    g = L[:, A].astype(np.int8) + L[:, B]
    H = 0.0
    mean = {}
    for code, name in ((2, "pp"), (1, "cx"), (0, "ee")):
        m = g == code
        n = m.sum(1)
        H = H + (m * RANK).sum(1) ** 2 / n
        mean[name] = (m * V).sum(1) / n
    return H, mean


def main():
    species, A, B, V, obs = load()
    n = len(species)
    k = int(obs.sum())
    total = comb(n, k)
    RANK = stats.rankdata(V)

    H0, m0 = group_stats(obs[None, :], A, B, V, RANK)
    obs_H = float(H0[0])
    pp, cx, ee = m0["pp"][0], m0["cx"][0], m0["ee"][0]
    d_pp_ee, d_pp_cx, d_cx_ee = pp - ee, pp - cx, cx - ee

    print("%d species (%d pathogenic, %d environmental), %d pairs"
          % (n, k, n - k, len(V)))
    print("observed mean pn/ps    pathogenic-pathogenic %.4f   cross %.4f   "
          "environmental-environmental %.4f" % (pp, cx, ee))

    gobs = obs[A].astype(int) + obs[B]
    kw = stats.kruskal(V[gobs == 2], V[gobs == 1], V[gobs == 0])
    print("\nfor comparison, treating the 231 pairs as independent:")
    print("   Kruskal-Wallis H = %.2f, p = %.3g   <-- anti-conservative"
          % (kw.statistic, kw.pvalue))

    ge = c_pp_ee = c_pp_cx = c_cx_ee = seen = 0
    combos = itertools.combinations(range(n), k)
    while True:
        buf = np.zeros((CHUNK, n), dtype=bool)
        i = 0
        for combo in combos:
            buf[i, list(combo)] = True
            i += 1
            if i == CHUNK:
                break
        if i == 0:
            break
        H, m = group_stats(buf[:i], A, B, V, RANK)
        ge += int((H >= obs_H - 1e-9).sum())
        c_pp_ee += int(((m["pp"] - m["ee"]) <= d_pp_ee + 1e-12).sum())
        c_pp_cx += int(((m["pp"] - m["cx"]) <= d_pp_cx + 1e-12).sum())
        c_cx_ee += int(((m["cx"] - m["ee"]) <= d_cx_ee + 1e-12).sum())
        seen += i
        if i < CHUNK:
            break
    assert seen == total, "enumerated %d of %d relabellings" % (seen, total)

    out = [
        ("omnibus (three groups)", "", ge, ge / total),
        ("pathogenic-pathogenic vs environmental-environmental",
         d_pp_ee, c_pp_ee, c_pp_ee / total),
        ("pathogenic-pathogenic vs cross", d_pp_cx, c_pp_cx, c_pp_cx / total),
        ("cross vs environmental-environmental", d_cx_ee, c_cx_ee, c_cx_ee / total),
    ]

    print("\nexact permutation over all %d species-label assignments" % total)
    for label, diff, cnt, p in out:
        d = "" if diff == "" else "  (difference in means %+.4f)" % diff
        print("   %-52s %6d / %d   p = %.4g%s" % (label, cnt, total, p, d))

    RESULTS.mkdir(exist_ok=True)
    with open(RESULTS / "permutation_test.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["contrast", "difference_in_means", "n_as_extreme",
                    "n_assignments", "p_exact"])
        for label, diff, cnt, p in out:
            w.writerow([label, "" if diff == "" else "%.6f" % diff,
                        cnt, total, "%.6g" % p])
    print("\nwrote %s" % (RESULTS / "permutation_test.tsv"))


if __name__ == "__main__":
    main()
