#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""R3.3b — make the ascertainment-bias argument measurable instead of conceding
it blind.

One concern is that "no known human disease" may measure how much a
species has been studied rather than its biology. That is testable if study
effort is measured. Two objective, reproducible proxies per species, both from
NCBI:

  n_nuccore  number of nucleotide records for the organism
  n_pubmed   number of PubMed records matching the species name
  year_first earliest year attached to a nucleotide record for the organism

Then:
  1. Does study effort differ between the pathogenic and environmental sets?
     (If it does, that concern is confirmed as a real asymmetry in
     our own independent variable, and we should say so.)
  2. Does study effort predict the per-pair omega, and does the group contrast
     survive conditioning on it?
  3. The restriction route: restrict to COMPARABLY WELL-STUDIED species and
     re-run the contrast.

No value is hard-coded; everything written to results/ comes from the queries
logged in results/05_ncbi_raw.json.
"""
import json
import sys
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

labels = pd.read_csv(f"{A}/inputs/labels.tsv", sep="\t")
PATHOGENIC = set(labels.loc[labels.category == "pathogenic", "abbrev"])


def esearch_count(db, term):
    q = urllib.parse.urlencode({"db": db, "term": term, "retmode": "json",
                                "retmax": 0, "tool": TOOL})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as fh:
        d = json.load(fh)
    return int(d["esearchresult"]["count"])


def earliest_year(organism):
    """Oldest nucleotide record for the organism, by sorting on date."""
    q = urllib.parse.urlencode({
        "db": "nuccore", "term": f'"{organism}"[Organism]', "retmode": "json",
        "retmax": 1, "sort": "Date Released", "tool": TOOL})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as fh:
        d = json.load(fh)
    ids = d["esearchresult"].get("idlist") or []
    if not ids:
        return None
    q = urllib.parse.urlencode({"db": "nuccore", "id": ids[0],
                                "retmode": "json", "tool": TOOL})
    with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as fh:
        s = json.load(fh)
    rec = s["result"][ids[0]]
    for k in ("createdate", "updatedate"):
        v = rec.get(k)
        if v:
            return int(str(v)[:4])
    return None


raw = {}
rows = []
for r in labels.itertuples():
    sp, ab = r.species, r.abbrev
    try:
        nn = esearch_count("nuccore", f'"{sp}"[Organism]')
        time.sleep(0.4)
        np_ = esearch_count("pubmed", f'"{sp}"[All Fields]')
        time.sleep(0.4)
        # the older literature uses "<name> virus", not the ICTV binomial
        alt = sp.replace(" orthohantavirus", " virus")
        np_alt = esearch_count("pubmed", f'"{alt}"[All Fields]')
        time.sleep(0.4)
        yr = earliest_year(sp)
        time.sleep(0.4)
    except Exception as exc:                      # network trouble is possible
        print(f"QUERY FAILED for {ab}: {exc}", file=sys.stderr)
        nn = np_ = np_alt = yr = None
    raw[ab] = dict(species=sp, n_nuccore=nn, n_pubmed_binomial=np_,
                   n_pubmed_classic=np_alt, year_first_nuccore=yr)
    rows.append(dict(abbrev=ab, species=sp, category=r.category,
                     host_class=r.host_class, n_nuccore=nn,
                     n_pubmed=max(x for x in (np_, np_alt) if x is not None)
                     if (np_ is not None or np_alt is not None) else None,
                     year_first=yr))
    print(f"{ab:5s} nuccore={nn} pubmed={np_}/{np_alt} first={yr}")

with open(f"{A}/results/05_ncbi_raw.json", "w") as fh:
    json.dump({"queried_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "eutils": EUTILS, "records": raw}, fh, indent=1)

eff = pd.DataFrame(rows)
eff.to_csv(f"{A}/results/05_study_effort.tsv", sep="\t", index=False)

out = []
def say(s=""):
    out.append(s)
    print(s)

say("# R3.3b — study effort measured, not assumed")
say()
say(f"NCBI E-utilities queried {time.strftime('%Y-%m-%d', time.gmtime())}; raw")
say("counts in `results/05_ncbi_raw.json`, per-species table in")
say("`results/05_study_effort.tsv`. Counts are a snapshot and will drift, which")
say("is why the raw response is kept.")
say()
ok = eff.dropna(subset=["n_nuccore", "n_pubmed"])
say(f"Species with usable counts: {len(ok)} of {len(eff)}")
say()
say("| species | category | host | nuccore records | PubMed records | earliest nuccore year |")
say("|---|---|---|---|---|---|")
for r in eff.sort_values(["category", "abbrev"]).itertuples():
    say(f"| {r.abbrev} | {r.category} | {r.host_class} | {r.n_nuccore} | "
        f"{r.n_pubmed} | {r.year_first} |")
say()

if len(ok) >= 15:
    p = ok[ok.category == "pathogenic"]
    e = ok[ok.category == "environmental"]
    say("## 1. Is study effort asymmetric between the two groups?")
    say()
    say("| measure | pathogenic median | environmental median | Mann-Whitney U | p (two-sided) |")
    say("|---|---|---|---|---|")
    asym = {}
    for col in ("n_nuccore", "n_pubmed", "year_first"):
        a = p[col].dropna().values.astype(float)
        b = e[col].dropna().values.astype(float)
        if len(a) < 3 or len(b) < 3:
            continue
        U, pv = stats.mannwhitneyu(a, b, alternative="two-sided")
        asym[col] = dict(med_path=float(np.median(a)),
                         med_env=float(np.median(b)), U=float(U), p=float(pv))
        say(f"| {col} | {np.median(a):.0f} | {np.median(b):.0f} | {U:.1f} | "
            f"{pv:.4g} |")
    say()

    # 2/3. does the contrast survive restricting to comparably studied species?
    pair = pd.read_csv(f"{A}/inputs/pairwise_three_estimators.tsv", sep="\t")
    effmap = {r.abbrev: r for r in ok.itertuples()}

    def group(a, b):
        ga = a in PATHOGENIC
        gb = b in PATHOGENIC
        return "path-path" if (ga and gb) else (
            "env-env" if not (ga or gb) else "cross")

    pair["grp"] = [group(r.species_A, r.species_B) for r in pair.itertuples()]

    say("## 2-3. Restricting to comparably well-studied species")
    say()
    say("Thresholds are quantiles of the nuccore count across all 22 species, so")
    say("the restriction rule is stated in advance and not chosen to suit the")
    say("answer. Each row keeps only pairs in which BOTH species clear the")
    say("threshold.")
    say()
    say("| nuccore threshold | species kept (path/env) | pairs | path-path mean | cross mean | env-env mean | KW p | ordering preserved |")
    say("|---|---|---|---|---|---|---|---|")
    restr = {}
    counts = ok.set_index("abbrev")["n_nuccore"].astype(float)
    for q in (0.0, 0.25, 0.50):
        thr = float(counts.quantile(q))
        keepsp = set(counts[counts >= thr].index)
        npath = len(keepsp & PATHOGENIC)
        nenv = len(keepsp) - npath
        m = pair.apply(lambda r: r.species_A in keepsp and r.species_B in keepsp,
                       axis=1)
        sub = pair[m]
        arrs, means = [], {}
        for g in ("path-path", "cross", "env-env"):
            v = sub.loc[sub.grp == g, "yn_omega"].dropna().values
            arrs.append(v)
            means[g] = float(np.mean(v)) if len(v) else float("nan")
        if all(len(x) > 0 for x in arrs):
            H, pv = stats.kruskal(*arrs)
        else:
            H, pv = float("nan"), float("nan")
        order = means["path-path"] < means["cross"] < means["env-env"]
        restr[f"q{q}"] = dict(thr=thr, n_pairs=int(len(sub)), means=means,
                              p=float(pv), ordering=bool(order),
                              npath=npath, nenv=nenv)
        say(f"| >= {thr:.0f} (q{q:.2f}) | {len(keepsp)} ({npath}/{nenv}) | "
            f"{len(sub)} | {means['path-path']:.4f} | {means['cross']:.4f} | "
            f"{means['env-env']:.4f} | {pv:.3e} | {order} |")
    say()
    say("Estimator used above is YN00 omega, the best-corrected of the three.")
    say()
    with open(f"{A}/results/05_study_effort.json", "w") as fh:
        json.dump({"asymmetry": asym, "restricted": restr}, fh, indent=1)

with open(f"{A}/results/05_study_effort.md", "w") as fh:
    fh.write("\n".join(out) + "\n")
print("\nwrote results/05_study_effort.{md,tsv,json} and 05_ncbi_raw.json")
