#!/usr/bin/env python3
"""R2.3b step 4 - make the null result quantitative.

A non-significant p-value alone does not answer R2.3b. For each testable
species we report:
  - mean difference (human - animal) with a 95% bootstrap CI,
  - pooled within-group SD and Cohen's d,
  - the minimum difference that would have been detectable at 80% power,
    alpha 0.05, two-sided, at the realised group sizes (so a reader can see
    the comparison was not merely underpowered),
  - the same test on the manuscript's own CAI statistic, as a sensitivity check.
  - a de-duplicated (one sequence per country x year x group) re-test, because
    outbreak clusters make sequences non-independent.
"""
import collections
import csv
import json
import math
import os
import random
import sys

TSV = sys.argv[1]
OUTDIR = sys.argv[2]
MIN_N = 8
random.seed(20261012)

rows = list(csv.DictReader(open(TSV), delimiter="\t"))
for r in rows:
    r["cai_proper"] = float(r["cai_proper"])
    r["cai_manuscript_style"] = float(r["cai_manuscript_style"])


def mannwhitney_p(a, b):
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return None
    allv = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks = [0.0] * (n1 + n2)
    i = 0
    ties = []
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = r
        ties.append(j - i + 1)
        i = j + 1
    r1 = sum(ranks[k] for k in range(len(allv)) if allv[k][1] == 0)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u = min(u1, n1 * n2 - u1)
    n = n1 + n2
    tiesum = sum(t ** 3 - t for t in ties)
    sd = math.sqrt((n1 * n2 / 12.0) * ((n + 1) - tiesum / (n * (n - 1.0))))
    if sd == 0:
        return 1.0
    return math.erfc(abs((u - n1 * n2 / 2.0) / sd) / math.sqrt(2))


def mean(v):
    return sum(v) / len(v)


def sd(v):
    if len(v) < 2:
        return 0.0
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def boot_ci(a, b, reps=10000):
    d = []
    for _ in range(reps):
        d.append(mean([random.choice(a) for _ in a]) - mean([random.choice(b) for _ in b]))
    d.sort()
    return d[int(0.025 * reps)], d[int(0.975 * reps)]


by_sp = collections.defaultdict(lambda: collections.defaultdict(list))
for r in rows:
    by_sp[r["organism"]][r["group"]].append(r)

Z_A = 1.959963985  # two-sided alpha 0.05
Z_B = 0.8416212336  # 80% power

out = []
for sp, g in sorted(by_sp.items()):
    H, A = g.get("human", []), g.get("animal", [])
    if len(H) < MIN_N or len(A) < MIN_N:
        continue
    hp = [r["cai_proper"] for r in H]
    ap = [r["cai_proper"] for r in A]
    hm = [r["cai_manuscript_style"] for r in H]
    am = [r["cai_manuscript_style"] for r in A]
    sp_pooled = math.sqrt(((len(hp) - 1) * sd(hp) ** 2 + (len(ap) - 1) * sd(ap) ** 2) /
                          (len(hp) + len(ap) - 2))
    mde = (Z_A + Z_B) * sp_pooled * math.sqrt(1.0 / len(hp) + 1.0 / len(ap))
    lo, hi = boot_ci(hp, ap)

    # independence-trimmed re-test: one sequence per (country, year, group)
    def trim(recs):
        seen, keep = set(), []
        for r in recs:
            key = (r["country"].split(":")[0].strip(), (r["coll_date"] or "")[:4])
            if key in seen:
                continue
            seen.add(key)
            keep.append(r["cai_proper"])
        return keep
    th, ta = trim(H), trim(A)

    out.append({
        "organism": sp,
        "n_human": len(hp), "n_animal": len(ap),
        "mean_human": mean(hp), "mean_animal": mean(ap),
        "delta_mean": mean(hp) - mean(ap),
        "ci95_delta": [lo, hi],
        "sd_human": sd(hp), "sd_animal": sd(ap), "sd_pooled": sp_pooled,
        "cohens_d": (mean(hp) - mean(ap)) / sp_pooled if sp_pooled else None,
        "p_mwu_proper_cai": mannwhitney_p(hp, ap),
        "min_detectable_difference_80pct_power": mde,
        "p_mwu_manuscript_statistic": mannwhitney_p(hm, am),
        "mean_human_manuscript_statistic": mean(hm),
        "mean_animal_manuscript_statistic": mean(am),
        "trimmed_n_human": len(th), "trimmed_n_animal": len(ta),
        "p_mwu_trimmed": mannwhitney_p(th, ta) if len(th) >= 5 and len(ta) >= 5 else None,
    })

with open(os.path.join(OUTDIR, "R2.3b_power_sensitivity.json"), "w") as fh:
    json.dump(out, fh, indent=2, sort_keys=True)

print("%-32s %4s %4s %10s %22s %9s %9s %8s %8s" % (
    "organism", "nH", "nA", "delta", "95% CI delta", "d", "MDE80", "p_prop", "p_ms"))
for r in out:
    print("%-32s %4d %4d %+10.5f [%+9.5f,%+9.5f] %+9.3f %9.5f %8.3g %8.3g" % (
        r["organism"][:32], r["n_human"], r["n_animal"], r["delta_mean"],
        r["ci95_delta"][0], r["ci95_delta"][1], r["cohens_d"],
        r["min_detectable_difference_80pct_power"], r["p_mwu_proper_cai"],
        r["p_mwu_manuscript_statistic"]))
print()
for r in out:
    print("%-32s trimmed nH=%d nA=%d p=%s" % (
        r["organism"][:32], r["trimmed_n_human"], r["trimmed_n_animal"],
        ("%.3g" % r["p_mwu_trimmed"]) if r["p_mwu_trimmed"] is not None else "n/a"))
