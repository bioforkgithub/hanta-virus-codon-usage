#!/usr/bin/env python3
"""Phase 1: download all GenBank nuccore records for the six 'environmental-B'
orthohantavirus species, into one fresh directory per species.

Run with:  python3 -I 01_download_environmental_B.py <inputs_dir> <logs_dir>

Deterministic: no RNG is used in this phase; accessions are sorted before use.
Downloaded data is UNTRUSTED -- it is written under <inputs_dir>/<ABBREV>/ and
never executed or imported. This script lives outside that tree.
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = os.environ.get("ENTREZ_EMAIL", "you@example.org")
TOOL = "hantavirus_pinps_inventory"

# abbrev -> (organism query term, free-text cross-check term)
SPECIES = [
    ("OXBV", '"Oxbow orthohantavirus"[Organism]', '"Oxbow virus"[All Fields]'),
    ("RKPV", '"Rockport orthohantavirus"[Organism]', '"Rockport virus"[All Fields]'),
    ("SWSV", '"Seewis orthohantavirus"[Organism]', '"Seewis virus"[All Fields]'),
    ("CBNV", '"Cao Bang orthohantavirus"[Organism]', '"Cao Bang virus"[All Fields]'),
    ("QZNV", '"Quezon orthohantavirus"[Organism]', '"Quezon virus"[All Fields]'),
    ("JJUV", '"Jeju orthohantavirus"[Organism]', '"Jeju virus"[All Fields]'),
]

SLEEP = 0.40  # stay under the 3 req/s unauthenticated NCBI limit


def get(url, params, retries=4):
    params = dict(params)
    params.update({"email": EMAIL, "tool": TOOL})
    full = url + "?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(full, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as exc:  # transient NCBI 429/500
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"failed after {retries}: {full}: {last}")


def esearch(term):
    txt = get(f"{EUTILS}/esearch.fcgi",
              {"db": "nuccore", "term": term, "retmax": 5000, "retmode": "json"})
    js = json.loads(txt)["esearchresult"]
    return int(js["count"]), sorted(js.get("idlist", [])), js.get("querytranslation", "")


def efetch_gb(uids, out_path, log):
    """Fetch full GenBank flatfiles in batches, concatenated to out_path."""
    with open(out_path, "w", encoding="utf-8") as fh:
        for i in range(0, len(uids), 150):
            chunk = uids[i:i + 150]
            txt = get(f"{EUTILS}/efetch.fcgi",
                      {"db": "nuccore", "id": ",".join(chunk),
                       "rettype": "gbwithparts", "retmode": "text"})
            fh.write(txt)
            if not txt.endswith("\n"):
                fh.write("\n")
            log.append(f"    fetched {i + len(chunk)}/{len(uids)}")
            time.sleep(SLEEP)


def main():
    inputs_dir, logs_dir = sys.argv[1], sys.argv[2]
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    manifest, log = {}, []

    for abbrev, org_term, text_term in SPECIES:
        sp_dir = os.path.join(inputs_dir, abbrev)
        os.makedirs(sp_dir, exist_ok=True)
        os.chmod(sp_dir, 0o700)  # private by default (global rule)

        n_org, ids_org, qt_org = esearch(org_term)
        time.sleep(SLEEP)
        n_txt, ids_txt, qt_txt = esearch(text_term)
        time.sleep(SLEEP)

        extra = sorted(set(ids_txt) - set(ids_org))
        uids = sorted(set(ids_org) | set(ids_txt), key=int)
        log.append(f"{abbrev}: organism n={n_org} ({qt_org}); "
                   f"text n={n_txt} ({qt_txt}); extra-from-text={len(extra)}; "
                   f"union={len(uids)}")

        gb_path = os.path.join(sp_dir, "records.gb")
        efetch_gb(uids, gb_path, log)

        manifest[abbrev] = {
            "organism_term": org_term, "organism_query_translation": qt_org,
            "organism_count": n_org, "text_term": text_term,
            "text_query_translation": qt_txt, "text_count": n_txt,
            "uids_from_text_only": extra, "n_union_uids": len(uids),
            "gb_path": gb_path,
        }
        print(f"{abbrev}: union {len(uids)} records -> {gb_path}", flush=True)

    with open(os.path.join(inputs_dir, "download_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    with open(os.path.join(logs_dir, "01_download.log"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("done")


if __name__ == "__main__":
    main()
