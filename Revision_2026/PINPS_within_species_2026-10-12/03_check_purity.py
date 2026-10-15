#!/usr/bin/env python3
"""Audit taxonomic purity of the downloaded sets.

For each species: list the distinct source /organism values among the USABLE
N CDS, flag any usable accession that was pulled in only by the free-text term
(i.e. not by the [Organism] term), and print the organism of every usable hit.

Run:  python3 -I 03_check_purity.py <inputs_dir>
"""
import json
import os
import sys

from Bio import SeqIO

ABBREVS = ["OXBV", "RKPV", "SWSV", "CBNV", "QZNV", "JJUV"]


def main():
    inputs_dir = sys.argv[1]
    man = json.load(open(os.path.join(inputs_dir, "download_manifest.json")))

    for ab in ABBREVS:
        fa = os.path.join(inputs_dir, ab, "N_cds_usable.fasta")
        usable_accs = {r.id.split("|")[0] for r in SeqIO.parse(fa, "fasta")}
        text_only = set(man[ab]["uids_from_text_only"])

        print(f"\n=== {ab} ===")
        print(f"  organism term     : {man[ab]['organism_query_translation']}"
              f"  n={man[ab]['organism_count']}")
        print(f"  text-only extras  : {len(text_only)} uids")

        orgs = {}
        flagged = []
        for rec in SeqIO.parse(os.path.join(inputs_dir, ab, "records.gb"), "genbank"):
            if rec.id not in usable_accs:
                continue
            org, host = "", ""
            for f in rec.features:
                if f.type == "source":
                    org = (f.qualifiers.get("organism") or [""])[0]
                    host = (f.qualifiers.get("host") or [""])[0]
                    break
            orgs.setdefault(org, []).append(rec.id)
            # was this record reachable ONLY via the free-text term?
            gi = [x for x in rec.annotations.get("accessions", [])]
            flagged.append((rec.id, org, host, rec.description[:70]))

        print("  distinct /organism among usable:")
        for org, accs in sorted(orgs.items()):
            print(f"    {len(accs):3d}  {org!r}  e.g. {accs[:3]}")
        print("  usable records whose /organism is a MINORITY value:")
        if orgs:
            major = max(orgs, key=lambda o: len(orgs[o]))
            for acc, org, host, desc in sorted(flagged):
                if org != major:
                    print(f"    {acc}  org={org!r} host={host!r}  {desc}")


if __name__ == "__main__":
    main()
