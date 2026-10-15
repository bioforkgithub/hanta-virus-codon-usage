#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""R3.3b, part 3 — the vintage-matched restriction with GROUP SIZES shown.

Script 06 reported that the three-group ordering breaks when species are
restricted to comparable vintage. Before that can be stated as a result it
has to be clear how many pairs are left in each group, because if the
environmental-environmental cell collapses to one or two pairs the
restriction is uninformative rather than contradictory. Uses the years
already cached in results/06_first_year.json — no new NCBI queries.
"""
import json
import itertools

import numpy as np
import pandas as pd
from scipy import stats

A = ("/PATH/TO/"
     "Project_Hantavirus/Reviews_ans/analysis/R3.1_deconfound_2026-10-12")

years = json.load(open(f"{A}/results/06_first_year.json"))["years"]
eff = pd.read_csv(f"{A}/results/05_study_effort.tsv", sep="\t")
PATH = set(eff.loc[eff.category == "pathogenic", "abbrev"])
pair = pd.read_csv(f"{A}/inputs/pairwise_three_estimators.tsv", sep="\t")

def grp(a, b):
    ga, gb = a in PATH, b in PATH
    return "path-path" if (ga and gb) else ("env-env" if not (ga or gb)
                                            else "cross")
pair["grp"] = [grp(r.species_A, r.species_B) for r in pair.itertuples()]

out = []
def say(s=""):
    out.append(s)
    print(s)

say("# R3.3b, part 3 — vintage-matched restriction, with cell counts")
say()
say("Each row keeps only species whose earliest sequence record is at or before")
say("the cutoff, i.e. species that have had at least as long to be observed as")
say("each other. `n` columns are pair counts, which is the point: the")
say("restriction is only interpretable where the cells are populated.")
say()
for est, lab in (("yn_omega", "YN00 omega"),
                 ("ms_omega", "manuscript pn/ps")):
    say(f"## {lab}")
    say()
    say("| cutoff year | species (path/env) | n path-path | n cross | n env-env | "
        "mean path-path | mean cross | mean env-env | KW p | ordering preserved |")
    say("|---|---|---|---|---|---|---|---|---|---|")
    for cut in (1994, 1997, 1999, 2004, 2008, 2016):
        keep = {a for a, y in years.items() if y and y <= cut}
        npath = len(keep & PATH)
        nenv = len(keep) - npath
        m = pair.apply(lambda r: r.species_A in keep and r.species_B in keep,
                       axis=1)
        sub = pair[m]
        ns, means, arrs = {}, {}, []
        for g in ("path-path", "cross", "env-env"):
            v = sub.loc[sub.grp == g, est].dropna().values
            ns[g] = len(v)
            means[g] = float(np.mean(v)) if len(v) else float("nan")
            arrs.append(v)
        if all(len(x) > 1 for x in arrs):
            H, p = stats.kruskal(*arrs)
            pstr = f"{p:.3e}"
        else:
            pstr = "not testable"
        order = (means["path-path"] < means["cross"] < means["env-env"])
        say(f"| <= {cut} | {len(keep)} ({npath}/{nenv}) | {ns['path-path']} | "
            f"{ns['cross']} | {ns['env-env']} | {means['path-path']:.4f} | "
            f"{means['cross']:.4f} | {means['env-env']:.4f} | {pstr} | "
            f"{order} |")
    say()

say("## The three confounders are not independent of each other")
say()
old = sorted([a for a, y in years.items() if y and y <= 1994])
env_old = [a for a in old if a not in PATH]
say(f"Species with a sequence record by 1994: {', '.join(old)}.")
say(f"Of these, the non-pathogenic ones are: {', '.join(env_old)} — "
    f"{len(env_old)} species, and **both are rodent-borne**.")
say()
say("This is the core difficulty and it should be stated in the paper in these")
say("terms: the only environmental species that have been available long enough")
say("to be comparably well studied are the only two environmental species that")
say("are rodent-borne. Study effort, host class and disease status are")
say("therefore mutually confounded, not three separable covariates, and no")
say("restriction of this dataset can break all three at once. The")
say("environmental-environmental cell in any vintage-matched restriction is the")
say("single PHV-TULV pair, so the apparent reversal of the ordering in those")
say("rows rests on one pair and is NOT evidence that the effect reverses; it is")
say("evidence that the question cannot be answered with these 22 species.")
say()

with open(f"{A}/results/07_vintage_detail.md", "w") as fh:
    fh.write("\n".join(out) + "\n")
print("\nwrote results/07_vintage_detail.md")
