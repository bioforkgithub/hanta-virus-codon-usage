#!/usr/bin/env python3
"""R2.3b step 1 - harvest host metadata for orthohantavirus nucleotide records.

Feasibility first: we only want to know how many records carry an explicit
/host (or isolation_source) qualifier AND a usable protein-coding CDS, split
into human-derived vs animal/reservoir-derived.

No numbers are invented here: everything written to results/ comes from the
GenBank records downloaded in this run.
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = sys.argv[1]  # inputs dir for raw downloads
TOOL = {"tool": "hantavirus_review_R2.3b", "email": os.environ.get("ENTREZ_EMAIL", "you@example.org")}

# One taxon per query. Within-species comparison is the only design in which a
# human-vs-animal CAI contrast is not confounded by viral species identity.
TAXA = [
    "Andes orthohantavirus",
    "Sin Nombre orthohantavirus",
    "Puumala orthohantavirus",
    "Hantaan orthohantavirus",
    "Seoul orthohantavirus",
    "Dobrava-Belgrade orthohantavirus",
    "Tula orthohantavirus",
    "Prospect Hill orthohantavirus",
    "Bayou orthohantavirus",
    "Black Creek Canal orthohantavirus",
    "Choclo orthohantavirus",
    "Laguna Negra orthohantavirus",
]


def get(url):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return fh.read()
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write("retry %d: %s\n" % (attempt, exc))
            time.sleep(3 + 3 * attempt)
    raise SystemExit("failed: %s" % url)


def esearch(term):
    q = dict(TOOL)
    q.update(db="nuccore", term=term, retmax="5000", retmode="json")
    data = json.loads(get(BASE + "/esearch.fcgi?" + urllib.parse.urlencode(q)))
    return data["esearchresult"].get("idlist", []), data["esearchresult"]["count"]


def efetch(ids, path):
    with open(path, "wb") as out:
        for i in range(0, len(ids), 150):
            chunk = ids[i:i + 150]
            q = dict(TOOL)
            q.update(db="nuccore", id=",".join(chunk), rettype="gb", retmode="text")
            out.write(get(BASE + "/efetch.fcgi?" + urllib.parse.urlencode(q)))
            time.sleep(0.5)


manifest = {}
for taxon in TAXA:
    term = '"%s"[Organism]' % taxon
    ids, count = esearch(term)
    slug = taxon.split()[0].lower()
    path = os.path.join(OUT, "gb_%s.gb" % slug)
    manifest[taxon] = {"count_reported": count, "ids_retrieved": len(ids), "file": path}
    sys.stderr.write("%-34s count=%s retrieved=%d\n" % (taxon, count, len(ids)))
    if ids:
        efetch(ids, path)
    time.sleep(0.5)

with open(os.path.join(OUT, "harvest_manifest.json"), "w") as fh:
    json.dump(manifest, fh, indent=2)
print(json.dumps(manifest, indent=2))
