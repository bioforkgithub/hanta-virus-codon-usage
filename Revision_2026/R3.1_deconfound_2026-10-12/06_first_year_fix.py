#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""R3.3b, correction. In script 05 the `year_first` column is WRONG: NCBI's
`sort=Date Released` returns the most recent record first, so the value stored
was the newest record's year, not the oldest (HTNV came back as 2026). That
column must not be used, and is recomputed here.

Correct method: ask esearch how many records exist with publication date at or
before year Y, for each Y, and take the first Y with a non-zero count. One
count query per year per species, which is exact and needs no record download.

The resulting "year of earliest sequence record" is then used as the
second study-effort covariate (years since description) in the R3.3b analysis.
"""
import json
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

A = ("/PATH/TO/"
     "Project_Hantavirus/Reviews_ans/analysis/R3.1_deconfound_2026-10-12")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "hantavirus_review_R3.3b"


def count_upto(organism, year):
    q = urllib.parse.urlencode({
        "db": "nuccore", "term": f'"{organism}"[Organism]', "retmode": "json",
        "retmax": 0, "datetype": "pdat", "mindate": "1950",
        "maxdate": str(year), "tool": TOOL})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as fh:
        return int(json.load(fh)["esearchresult"]["count"])


def earliest(organism, lo=1975, hi=2026):
    """Binary search for the first year with any record."""
    if count_upto(organism, hi) == 0:
        return None
    time.sleep(0.4)
    if count_upto(organism, lo) > 0:
        return lo
    while hi - lo > 1:
        mid = (lo + hi) // 2
        time.sleep(0.4)
        if count_upto(organism, mid) > 0:
            hi = mid
        else:
            lo = mid
    return hi


eff = pd.read_csv(f"{A}/results/05_study_effort.tsv", sep="\t")
years = {}
for r in eff.itertuples():
    try:
        y = earliest(r.species)
    except Exception as exc:
        print(f"FAILED {r.abbrev}: {exc}")
        y = None
    years[r.abbrev] = y
    print(f"{r.abbrev:5s} earliest sequence-record year = {y}")
    time.sleep(0.4)

eff["year_first"] = [years[a] for a in eff.abbrev]        # overwrite the bad one
eff["years_since_first"] = [2026 - years[a] if years[a] else None
                            for a in eff.abbrev]
eff.to_csv(f"{A}/results/05_study_effort.tsv", sep="\t", index=False)
with open(f"{A}/results/06_first_year.json", "w") as fh:
    json.dump({"method": "esearch pdat count<=Y, binary search",
               "queried_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "years": years}, fh, indent=1)

out = []
def say(s=""):
    out.append(s)
    print(s)

say("# R3.3b, part 2 — years since first sequence record")
say()
say("`year_first` as written by script 05 was wrong (NCBI's `sort=Date")
say("Released` is newest-first), and has been recomputed here by exact")
say("date-bounded counts. Raw values in `results/06_first_year.json`.")
say()
say("This is the year of the earliest *sequence record*, which is a lower bound")
say("on how long the species has been available to study, not the year of")
say("taxonomic description. We say so in the manuscript rather than calling it")
say("a description date.")
say()
ok = eff.dropna(subset=["year_first"])
say("| species | category | host | earliest record | years since | nuccore | PubMed |")
say("|---|---|---|---|---|---|---|")
for r in ok.sort_values(["category", "year_first"]).itertuples():
    say(f"| {r.abbrev} | {r.category} | {r.host_class} | {int(r.year_first)} | "
        f"{int(r.years_since_first)} | {r.n_nuccore} | {r.n_pubmed} |")
say()
p = ok[ok.category == "pathogenic"]["year_first"].astype(float).values
e = ok[ok.category == "environmental"]["year_first"].astype(float).values
U, pv = stats.mannwhitneyu(p, e, alternative="two-sided")
say(f"Earliest-record year, pathogenic median {np.median(p):.0f} "
    f"(range {p.min():.0f}-{p.max():.0f}) vs environmental median "
    f"{np.median(e):.0f} (range {e.min():.0f}-{e.max():.0f}); "
    f"Mann-Whitney U = {U:.1f}, p = {pv:.4g}.")
say()

# does the group contrast survive restricting to species of similar vintage?
pair = pd.read_csv(f"{A}/inputs/pairwise_three_estimators.tsv", sep="\t")
PATH = set(eff.loc[eff.category == "pathogenic", "abbrev"])
def grp(a, b):
    ga, gb = a in PATH, b in PATH
    return "path-path" if (ga and gb) else ("env-env" if not (ga or gb)
                                            else "cross")
pair["grp"] = [grp(r.species_A, r.species_B) for r in pair.itertuples()]
yrs = ok.set_index("abbrev")["year_first"].astype(float)

say("## Restricting to species of comparable vintage")
say()
say("| cutoff | species kept (path/env) | pairs | path-path | cross | env-env | KW p | ordering preserved |")
say("|---|---|---|---|---|---|---|---|")
res = {}
for cut in sorted({int(yrs.quantile(q)) for q in (0.0, 0.25, 0.5)}):
    keep = set(yrs[yrs <= cut].index) if cut > yrs.min() else set(yrs.index)
    keep = set(yrs[yrs <= cut].index)
    npath = len(keep & PATH)
    nenv = len(keep) - npath
    m = pair.apply(lambda r: r.species_A in keep and r.species_B in keep, axis=1)
    sub = pair[m]
    arrs, means = [], {}
    for g in ("path-path", "cross", "env-env"):
        v = sub.loc[sub.grp == g, "yn_omega"].dropna().values
        arrs.append(v)
        means[g] = float(np.mean(v)) if len(v) else float("nan")
    if all(len(x) > 1 for x in arrs):
        H, pv2 = stats.kruskal(*arrs)
    else:
        H, pv2 = float("nan"), float("nan")
    order = means["path-path"] < means["cross"] < means["env-env"]
    res[cut] = dict(n_pairs=int(len(sub)), means=means, p=float(pv2),
                    ordering=bool(order), npath=npath, nenv=nenv)
    say(f"| earliest record <= {cut} | {len(keep)} ({npath}/{nenv}) | "
        f"{len(sub)} | {means['path-path']:.4f} | {means['cross']:.4f} | "
        f"{means['env-env']:.4f} | {pv2:.3e} | {order} |")
say()
say("Estimator: YN00 omega.")
say()
with open(f"{A}/results/06_vintage_restricted.json", "w") as fh:
    json.dump(res, fh, indent=1, default=str)
with open(f"{A}/results/06_first_year.md", "w") as fh:
    fh.write("\n".join(out) + "\n")
print("\nwrote results/06_first_year.{md,json}, 06_vintage_restricted.json; "
      "05_study_effort.tsv year_first column corrected")
