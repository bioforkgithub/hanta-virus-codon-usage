#!/usr/bin/env python3
"""R2.3b step 2 - parse GenBank, classify isolate host, extract complete N CDS.

Inclusion rules (all applied in code, all exclusions counted):
  1. record carries an explicit /host or /isolation_source qualifier that can be
     resolved to HUMAN or to a non-human ANIMAL species;
  2. record carries a CDS annotated as nucleocapsid / nucleoprotein / N protein
     (the manuscript's S-segment gene);
  3. that CDS is complete: length %% 3 == 0, >= 1200 nt, no /pseudo, not
     partial (no '<' or '>' in the location), ACGT only;
  4. one sequence per isolate/strain per taxon (first occurrence kept), so that
     resequenced strains do not inflate n.
"""
import csv
import glob
import json
import os
import re
import sys

from Bio import SeqIO

RAW = sys.argv[1]
OUTDIR = sys.argv[2]

HUMAN_PAT = re.compile(
    r"homo sapiens|^human$|human\b|patient|\bhfrs\b|hantavirus pulmonary syndrome"
    r"|hantavirus cardiopulmonary|\bhcps\b|\bhps\b|nephropathia",
    re.I,
)
# Reservoir / non-human host genera seen in orthohantavirus records.
ANIMAL_GENERA = {
    "oligoryzomys", "peromyscus", "apodemus", "rattus", "myodes", "clethrionomys",
    "microtus", "sigmodon", "oryzomys", "calomys", "akodon", "abrothrix",
    "necromys", "loxodontomys", "eligmodontia", "suncus", "crocidura", "sorex",
    "neomys", "talpa", "urotrichus", "neurotrichus", "scalopus", "blarina",
    "anourosorex", "episoriculus", "chaerephon", "rousettus", "mus",
    "bandicota", "niviventer", "apodemus/", "eothenomys", "lemmus", "myopus",
    "dicrostonyx", "ondatra", "arvicola", "reithrodontomys", "zygodontomys",
    "holochilus", "nectomys", "bolomys", "graomys", "phyllotis", "chelemys",
    "hylaeamys", "cerradomys", "sooretamys", "euryoryzomys", "oecomys",
    "neacomys", "rhipidomys", "thaptomys", "juliomys", "cricetulus",
    "erinaceus", "soricidae", "neodon", "lasiopodomys", "alexandromys",
}
N_PAT = re.compile(r"nucleocapsid|nucleoprotein|\bN protein\b|^N$|^NP$", re.I)


def classify(host_strings):
    """Return ('human'|'animal'|None, evidence string)."""
    for s in host_strings:
        s = s.strip()
        if not s:
            continue
        low = s.lower()
        g = re.split(r"[ _/,;]", low)[0]
        if g in ANIMAL_GENERA:
            return "animal", s
        if HUMAN_PAT.search(low):
            return "human", s
    for s in host_strings:
        low = s.lower()
        if any(gen in low for gen in ANIMAL_GENERA):
            return "animal", s
    return None, ";".join(x for x in host_strings if x)


rows = []
excl = {
    "no_host_qualifier": 0,
    "host_unresolved": 0,
    "no_N_cds": 0,
    "N_cds_partial_or_bad": 0,
    "duplicate_isolate": 0,
}
host_vocab = {}
seen = set()
per_taxon = {}

for path in sorted(glob.glob(os.path.join(RAW, "gb_*.gb"))):
    for rec in SeqIO.parse(path, "genbank"):
        src = [f for f in rec.features if f.type == "source"]
        hosts, isol, strain, segment, organism = [], [], "", "", ""
        organism = rec.annotations.get("organism", "")
        if src:
            q = src[0].qualifiers
            hosts = q.get("host", [])
            isol = q.get("isolation_source", [])
            strain = (q.get("strain") or q.get("isolate") or [""])[0]
            segment = (q.get("segment") or [""])[0]
        per_taxon.setdefault(organism, {"records": 0, "with_host": 0, "human": 0, "animal": 0, "usable": 0})
        per_taxon[organism]["records"] += 1
        for h in hosts + isol:
            host_vocab[h] = host_vocab.get(h, 0) + 1
        if not (hosts or isol):
            excl["no_host_qualifier"] += 1
            continue
        per_taxon[organism]["with_host"] += 1
        grp, evid = classify(hosts + isol)
        if grp is None:
            excl["host_unresolved"] += 1
            continue
        per_taxon[organism][grp] += 1

        cds = None
        for f in rec.features:
            if f.type != "CDS":
                continue
            names = f.qualifiers.get("product", []) + f.qualifiers.get("gene", []) + \
                f.qualifiers.get("note", [])
            if any(N_PAT.search(n) for n in names):
                cds = f
                break
        if cds is None:
            excl["no_N_cds"] += 1
            continue
        loc = str(cds.location)
        if "<" in loc or ">" in loc or "pseudo" in cds.qualifiers:
            excl["N_cds_partial_or_bad"] += 1
            continue
        try:
            seq = str(cds.extract(rec.seq)).upper()
        except Exception:  # noqa: BLE001
            excl["N_cds_partial_or_bad"] += 1
            continue
        if len(seq) % 3 or len(seq) < 1200 or set(seq) - set("ACGT"):
            excl["N_cds_partial_or_bad"] += 1
            continue

        key = (organism, strain.lower() or rec.id)
        if key in seen:
            excl["duplicate_isolate"] += 1
            continue
        seen.add(key)
        per_taxon[organism]["usable"] += 1
        rows.append({
            "accession": rec.id,
            "organism": organism,
            "segment": segment,
            "strain": strain,
            "group": grp,
            "host_evidence": evid,
            "country": (src[0].qualifiers.get("country") or
                        src[0].qualifiers.get("geo_loc_name") or [""])[0] if src else "",
            "coll_date": (src[0].qualifiers.get("collection_date") or [""])[0] if src else "",
            "cds_len": len(seq),
            "seq": seq,
        })

with open(os.path.join(OUTDIR, "R2.3b_dataset.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=[k for k in rows[0] if k != "seq"], delimiter="\t") \
        if rows else None
    if w:
        w.writeheader()
        for r in rows:
            w.writerow({k: v for k, v in r.items() if k != "seq"})
with open(os.path.join(OUTDIR, "R2.3b_cds.fasta"), "w") as fh:
    for r in rows:
        fh.write(">%s|%s|%s|%s\n%s\n" % (r["accession"], r["organism"].replace(" ", "_"),
                                         r["group"], r["strain"].replace(" ", "_"), r["seq"]))

summary = {"exclusions": excl, "per_taxon": per_taxon, "n_usable": len(rows)}
with open(os.path.join(OUTDIR, "R2.3b_parse_summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2, sort_keys=True)
with open(os.path.join(OUTDIR, "host_vocabulary.tsv"), "w") as fh:
    for h, c in sorted(host_vocab.items(), key=lambda kv: -kv[1]):
        fh.write("%d\t%s\n" % (c, h))

print(json.dumps({"n_usable": len(rows), "exclusions": excl}, indent=2))
for org, d in sorted(per_taxon.items()):
    print("%-36s recs=%5d host=%4d human=%4d animal=%4d usable=%4d" %
          (org[:36], d["records"], d["with_host"], d["human"], d["animal"], d["usable"]))
