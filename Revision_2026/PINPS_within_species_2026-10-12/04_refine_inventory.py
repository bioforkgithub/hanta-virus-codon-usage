#!/usr/bin/env python3
"""Phase 3: correct three defects found by the purity audit, and re-inventory.

  (a) free-text query contamination -- usable hits whose /organism belongs to a
      DIFFERENT species (e.g. Cao Bang records pulled into the SWSV set);
  (b) RefSeq redundancy -- NC_* records are curated copies of an original
      submission, so an NC_/original pair is ONE observation, not two;
  (c) divergent sub-lineages inside one NCBI species taxon (e.g. Lianghe virus
      inside Orthohantavirus caobangense), which would turn 'within-species
      polymorphism' back into between-lineage divergence.

Also reports within-species pairwise divergence, which is the quantity that has
to be unsaturated for piN/piS to be estimable at all.

Run: python3 -I 04_refine_inventory.py <inputs_dir> <results_dir> <logs_dir>
Deterministic: no RNG (accessions sorted); clustering is single-linkage at a
fixed 0.95 identity threshold.
"""
import itertools
import json
import os
import re
import sys
from collections import Counter

from Bio import SeqIO

ABBREVS = ["OXBV", "RKPV", "SWSV", "CBNV", "QZNV", "JJUV"]
HOSTCLASS = {"OXBV": "mole", "RKPV": "shrew", "SWSV": "shrew",
             "CBNV": "shrew", "QZNV": "bat", "JJUV": "shrew"}
# the organism tokens that legitimately denote each target species
ACCEPT = {
    "OXBV": ["oxbow"],
    "RKPV": ["rockport"],
    "SWSV": ["seewis"],
    "CBNV": ["caobangense", "cao bang", "lianghe"],
    "QZNV": ["quezon"],
    "JJUV": ["jeju", "jejuense"],
}
CLUSTER_ID = 0.95   # single-linkage nt-identity threshold for sub-lineages
MIN_USABLE = 10
MIN_HAPLOTYPES = 4


def pdist(a, b):
    """Ungapped p-distance over sites where both are unambiguous ACGT."""
    n = comp = 0
    for x, y in zip(a, b):
        if x in "ACGT" and y in "ACGT":
            comp += 1
            if x != y:
                n += 1
    return (n / comp if comp else None), comp


def clusters(accs, seqs, thresh):
    """Single-linkage clusters at >= thresh identity (equal-length pairs only)."""
    parent = {a: a for a in accs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in itertools.combinations(accs, 2):
        if len(seqs[a]) != len(seqs[b]):
            continue
        d, comp = pdist(seqs[a], seqs[b])
        if d is not None and (1 - d) >= thresh:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    groups = {}
    for a in accs:
        groups.setdefault(find(a), []).append(a)
    return sorted((sorted(v) for v in groups.values()), key=len, reverse=True)


def main():
    inputs_dir, results_dir, logs_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    report, log = {}, []

    for ab in ABBREVS:
        # ---- load usable CDS and their records ----
        fa = os.path.join(inputs_dir, ab, "N_cds_usable.fasta")
        usable = {r.id.split("|")[0]: str(r.seq).upper()
                  for r in SeqIO.parse(fa, "fasta")}
        meta = {}
        for rec in SeqIO.parse(os.path.join(inputs_dir, ab, "records.gb"), "genbank"):
            if rec.id not in usable:
                continue
            q = {}
            for f in rec.features:
                if f.type == "source":
                    q = f.qualifiers
                    break
            meta[rec.id] = {
                "organism": (q.get("organism") or [""])[0],
                "host": (q.get("host") or [""])[0],
                "country": (q.get("country") or q.get("geo_loc_name") or [""])[0],
                "isolation_source": (q.get("isolation_source") or [""])[0],
                "strain": (q.get("strain") or q.get("isolate") or [""])[0],
                "desc": rec.description,
                "is_refseq": rec.id.startswith("NC_"),
            }

        # ---- (a) drop wrong-species records ----
        wrong = sorted(a for a in usable
                       if not any(t in meta[a]["organism"].lower()
                                  for t in ACCEPT[ab]))
        for a in wrong:
            usable.pop(a)

        # ---- (b) collapse RefSeq/original identical pairs ----
        refseq_dropped = []
        for a in sorted(a for a in usable if meta[a]["is_refseq"]):
            twin = sorted(b for b, s in usable.items()
                          if b != a and s == usable[a] and not meta[b]["is_refseq"])
            if twin:
                refseq_dropped.append({"refseq": a, "original_kept": twin[0]})
                usable.pop(a)

        accs = sorted(usable)
        n_ind = len(accs)

        # ---- (c) sub-lineage structure + divergence ----
        cl = clusters(accs, usable, CLUSTER_ID) if n_ind >= 2 else [accs]
        ds, skipped = [], 0
        for a, b in itertools.combinations(accs, 2):
            if len(usable[a]) != len(usable[b]):
                skipped += 1
                continue
            d, comp = pdist(usable[a], usable[b])
            if d is not None:
                ds.append(d)
        ds.sort()

        # ---- duplicates among the independent set ----
        hap = Counter(usable[a] for a in accs)
        n_hap = len(hap)

        orgs = Counter(meta[a]["organism"] for a in accs)
        hosts = sorted({meta[a]["host"] for a in accs if meta[a]["host"]})
        ctys = sorted({meta[a]["country"].split(":")[0].strip()
                       for a in accs if meta[a]["country"]})
        nh = sum(1 for a in accs if meta[a]["host"])
        nc = sum(1 for a in accs if meta[a]["country"])
        ni = sum(1 for a in accs if meta[a]["isolation_source"])

        # largest single sub-lineage = the defensible within-lineage sample
        n_largest = len(cl[0]) if cl else 0
        sufficient = (n_ind >= MIN_USABLE and n_hap >= MIN_HAPLOTYPES)
        sufficient_strict = (n_largest >= MIN_USABLE and
                             len(Counter(usable[a] for a in cl[0])) >= MIN_HAPLOTYPES
                             if cl else False)

        report[ab] = {
            "host_class": HOSTCLASS[ab],
            "n_usable_before_refinement": len(accs) + len(wrong) + len(refseq_dropped),
            "wrong_species_removed": [{"acc": a, "organism": meta[a]["organism"]}
                                      for a in wrong],
            "refseq_redundancy_removed": refseq_dropped,
            "n_independent": n_ind,
            "accessions": accs,
            "n_distinct_haplotypes": n_hap,
            "n_duplicate_sequences": n_ind - n_hap,
            "organism_labels": dict(orgs),
            "sublineage_clusters_at_95pct_id": [len(c) for c in cl],
            "sublineage_members": {f"cluster{i+1}": c for i, c in enumerate(cl)},
            "n_largest_sublineage": n_largest,
            "pairwise_pdist_min": round(ds[0], 5) if ds else None,
            "pairwise_pdist_median": round(ds[len(ds) // 2], 5) if ds else None,
            "pairwise_pdist_max": round(ds[-1], 5) if ds else None,
            "n_pairs_compared": len(ds), "n_pairs_skipped_length": skipped,
            "frac_with_host": round(nh / n_ind, 3) if n_ind else None,
            "frac_with_country": round(nc / n_ind, 3) if n_ind else None,
            "frac_with_isolation_source": round(ni / n_ind, 3) if n_ind else None,
            "hosts": hosts, "countries": ctys, "n_countries": len(ctys),
            "sufficient_whole_taxon": sufficient,
            "sufficient_largest_sublineage": bool(sufficient_strict),
        }
        log.append(f"{ab}: usable_raw->{report[ab]['n_usable_before_refinement']} "
                   f"wrong_species={len(wrong)} refseq_dropped={len(refseq_dropped)} "
                   f"independent={n_ind} haplotypes={n_hap} "
                   f"clusters={[len(c) for c in cl]} "
                   f"pdist_max={report[ab]['pairwise_pdist_max']} "
                   f"sufficient={sufficient}")
        print(log[-1], flush=True)

        # refined FASTA
        with open(os.path.join(inputs_dir, ab, "N_cds_refined.fasta"), "w") as fh:
            for a in accs:
                m = meta[a]
                fh.write(f">{a}|{ab}|{m['strain'] or 'NA'}|host={m['host'] or 'NA'}"
                         f"|country={m['country'] or 'NA'}|org={m['organism']}\n"
                         f"{usable[a]}\n")

    with open(os.path.join(results_dir, "refined_inventory_environmental-B.json"), "w") as fh:
        json.dump({"cluster_identity_threshold": CLUSTER_ID,
                   "MIN_USABLE": MIN_USABLE, "MIN_HAPLOTYPES": MIN_HAPLOTYPES,
                   "species": report}, fh, indent=2, sort_keys=True)
    with open(os.path.join(logs_dir, "04_refine.log"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("ok")


if __name__ == "__main__":
    main()
