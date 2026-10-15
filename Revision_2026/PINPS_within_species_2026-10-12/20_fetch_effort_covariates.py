#!/usr/bin/env python3
import os  # noqa: E402  (added on publication, for ENTREZ_EMAIL)
"""Fetch study-effort covariates for the confound-direct statistical treatment.

For each species, find the EARLIEST YEAR in which any GenBank nuccore record
matching the phase-1/2 download query was published (datetype=pdat), by exact
binary search over maxdate.  This is a *measured proxy* for "years since species
description": the year the species first entered GenBank.  It is NOT the formal
taxonomic description year, which is not present in any file on disk; see the
results README.

Deterministic: pure binary search, no randomness.  The only non-determinism is
NCBI's own database contents, so the retrieval date is recorded in the output.

Usage: python3 -I 20_fetch_effort_covariates.py <ANALYSIS_ROOT> <LABELS_TSV> <OUT_JSON>
"""
import json
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
TERM = '("{species}"[Organism]) AND 300:5000[SLEN] NOT patent[Properties]'
SLEEP = 0.40  # stay under the 3 req/s anonymous limit


def esearch_count(species, mindate=None, maxdate=None, log=None):
    params = {
        "db": "nuccore",
        "term": TERM.format(species=species),
        "retmax": "0",
        "tool": "hantavirus_pinps_stats",
        "email": os.environ.get("ENTREZ_EMAIL", "you@example.org"),
    }
    if mindate is not None:
        params["datetype"] = "pdat"
        params["mindate"] = str(mindate)
        params["maxdate"] = str(maxdate)
    url = EUTILS + "?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                xml = fh.read()
            break
        except Exception as exc:  # noqa: BLE001
            if attempt == 3:
                raise
            if log:
                log.write("retry %d for %s (%s)\n" % (attempt, species, exc))
            time.sleep(2.0 * (attempt + 1))
    cnt = int(ET.fromstring(xml).findtext("Count"))
    if log:
        log.write("%-6s mindate=%s maxdate=%s -> %d\n"
                  % (species, mindate, maxdate, cnt))
        log.flush()
    time.sleep(SLEEP)
    return cnt


def earliest_year(species, lo=1975, hi=2026, log=None):
    """Smallest Y in [lo, hi] with >=1 record published in [lo, Y]."""
    if esearch_count(species, lo, hi, log) == 0:
        return None, 0
    total = esearch_count(species, lo, hi, log)
    a, b = lo, hi
    while a < b:
        mid = (a + b) // 2
        if esearch_count(species, lo, mid, log) > 0:
            b = mid
        else:
            a = mid + 1
    return a, total


def main():
    root, labels_tsv, out_json = sys.argv[1], sys.argv[2], sys.argv[3]
    want = set()
    rows = {}
    with open(labels_tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            rows[f["abbrev"]] = f
            want.add(f["abbrev"])

    out = {"retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "query_template": TERM,
           "search_window": [1975, 2026],
           "species": {}}
    with open(root + "/logs/20_fetch_effort_covariates.log", "w") as log:
        for ab in sorted(want):
            sp = rows[ab]["species"]
            yr, total = earliest_year(sp, log=log)
            out["species"][ab] = {"species": sp,
                                  "first_genbank_year": yr,
                                  "esearch_count_1975_2026": total}
            print(ab, sp, yr, total, flush=True)
    with open(out_json, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)


if __name__ == "__main__":
    main()
