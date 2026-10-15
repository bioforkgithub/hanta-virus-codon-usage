#!/usr/bin/env python3
"""Phase 5: assemble the final inventory TSV (+ ambiguity stats) for the six
'environmental-B' species.

Run: python3 -I 06_write_inventory.py <inputs_dir> <results_dir>
Deterministic: no RNG anywhere in the inventory; thresholds are constants.
"""
import json
import os
import re
import sys

from Bio import SeqIO

ABBREVS = ["OXBV", "RKPV", "SWSV", "CBNV", "QZNV", "JJUV"]
# sufficiency rule (see notes file for the justification)
MIN_INDEPENDENT = 10      # sequences
MIN_HAPLOTYPES = 4        # distinct haplotypes
MAX_LINEAGE_PDIST = 0.05  # a sample spanning more than this is not a population


def main():
    inputs_dir, results_dir = sys.argv[1], sys.argv[2]
    inv = json.load(open(os.path.join(
        results_dir, "inventory_environmental-B.json")))["species"]
    ref = json.load(open(os.path.join(
        results_dir, "refined_inventory_environmental-B.json")))["species"]
    st = json.load(open(os.path.join(
        results_dir, "structure_sensitivity_environmental-B.json")))["species"]

    rows = []
    for ab in ABBREVS:
        i, r, s = inv[ab], ref[ab], st[ab]
        b95 = s["by_threshold"]["0.95"]
        b90 = s["by_threshold"]["0.90"]

        # ambiguity over the refined set
        amb = []
        for rec in SeqIO.parse(os.path.join(inputs_dir, ab, "N_cds_refined.fasta"),
                               "fasta"):
            seq = str(rec.seq).upper()
            runs = re.findall(r"N+", seq)
            amb.append((sum(1 for c in seq if c not in "ACGT"),
                        max((len(x) for x in runs), default=0)))
        max_nonacgt = max((a for a, _ in amb), default=0)
        max_nrun = max((b for _, b in amb), default=0)

        pmax = r["pairwise_pdist_max"]
        count_ok = (r["n_independent"] >= MIN_INDEPENDENT and
                    r["n_distinct_haplotypes"] >= MIN_HAPLOTYPES)
        lineage_ok = (pmax is not None and pmax <= MAX_LINEAGE_PDIST)
        rows.append({
            "abbrev": ab, "category": "environmental",
            "host_class": r["host_class"],
            "n_genbank_records": i["n_genbank_records"],
            "n_complete_N_cds_raw": i["n_complete_N_cds"],
            "n_wrong_species_removed": len(r["wrong_species_removed"]),
            "n_refseq_redundancy_removed": len(r["refseq_redundancy_removed"]),
            "n_usable_independent": r["n_independent"],
            "n_distinct_haplotypes": r["n_distinct_haplotypes"],
            "n_duplicate_sequences": r["n_duplicate_sequences"],
            "cds_len_min": i["cds_len_min"], "cds_len_median": i["cds_len_median"],
            "cds_len_max": i["cds_len_max"],
            "max_nonACGT_in_usable": max_nonacgt,
            "max_N_run_in_usable": max_nrun,
            "n_excluded_no_N_cds": i["exclusions"].get("no_N_CDS_in_record", 0),
            "n_excluded_partial_cds": i["exclusions"].get("partial CDS", 0),
            "n_excluded_length": (i["exclusions"].get("too short", 0) +
                                  i["exclusions"].get("near-complete only", 0) +
                                  i["exclusions"].get("over-long", 0)),
            "n_excluded_ambiguity": 0,
            "pdist_median": r["pairwise_pdist_median"],
            "pdist_max": pmax,
            "n_lineages_at_95pct": b95["n_clusters"],
            "largest_lineage_n_at_95pct": b95["largest_n"],
            "largest_lineage_haplotypes_at_95pct": b95["largest_n_haplotypes"],
            "largest_lineage_pdist_max_at_95pct": b95["largest_pdist_max"],
            "n_lineages_at_90pct": b90["n_clusters"],
            "largest_lineage_n_at_90pct": b90["largest_n"],
            "largest_lineage_pdist_max_at_90pct": b90["largest_pdist_max"],
            "frac_with_host": r["frac_with_host"],
            "frac_with_country": r["frac_with_country"],
            "frac_with_isolation_source": r["frac_with_isolation_source"],
            "n_countries": r["n_countries"],
            "countries": ";".join(r["countries"]) or "NA",
            "hosts": ";".join(r["hosts"]) or "NA",
            "organism_labels": ";".join(f"{k}={v}" for k, v in
                                        sorted(r["organism_labels"].items())),
            "count_criterion_met": count_ok,
            "single_lineage_criterion_met": lineage_ok,
            "sufficient_for_pinps": count_ok,
            "sufficient_single_lineage": count_ok and lineage_ok,
        })

    cols = list(rows[0].keys())
    out = os.path.join(results_dir, "inventory_environmental-B.tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("wrote", out)
    for r in rows:
        print(f"  {r['abbrev']:5s} rec={r['n_genbank_records']:4d} "
              f"complete={r['n_complete_N_cds_raw']:3d} "
              f"indep={r['n_usable_independent']:3d} "
              f"hap={r['n_distinct_haplotypes']:3d} "
              f"pdistmax={r['pdist_max']} "
              f"lin95={r['n_lineages_at_95pct']} "
              f"largest95={r['largest_lineage_n_at_95pct']} "
              f"count_ok={r['count_criterion_met']} "
              f"lineage_ok={r['single_lineage_criterion_met']}")


if __name__ == "__main__":
    main()
