#!/usr/bin/env python3
"""
Assemble the per-species inventory TSV from the per-species summary.json +
exclusion_breakdown.json produced by fetch_N_cds.py / diagnose_excluded.py.

Usage: python3 -I build_inventory.py <INPUTS_DIR> <LABELS_TSV> <OUT_TSV> ABBREV [ABBREV...]

No RNG. Output row order = the ABBREV order given on the command line.
"""
import sys, os, json

# ---- sufficiency thresholds (stated and justified in the notes file) -------
# piN/piS at the gene level needs enough sequences that pS is non-zero and its
# sampling variance is not dominated by the pair count. With a 429-codon gene
# and within-species pS typically 0.02-0.06, n=20 gives 190 pairwise
# comparisons and a pS standard error of roughly 15-20% of the estimate --
# enough for a per-species point value in a 10-vs-12 group comparison. Below
# n=10 the per-species estimate is dominated by which few isolates happened to
# be sequenced, and a single divergent isolate can move piN/piS by a factor of
# two, so n<10 is reported as insufficient.
N_SUFFICIENT = 20
N_MARGINAL = 10

COLS = ["abbrev", "species", "category", "host_class", "reference_accession",
        "esearch_records", "records_M_or_L", "records_partial_N",
        "records_no_cds_feature", "records_unclassifiable",
        "excl_incomplete_cds", "excl_ambiguity",
        "n_complete_n_cds", "n_usable", "n_unique_usable",
        "n_duplicate_records", "largest_duplicate_group",
        "cds_len_min", "cds_len_modal", "cds_len_max",
        "host_available", "country_available", "date_available",
        "n_countries", "countries", "sufficient", "notes"]


def main():
    inputs, labels_tsv, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
    abbrevs = sys.argv[4:]

    lab = {}
    with open(labels_tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            lab[p[0]] = dict(zip(hdr, p))

    rows = []
    for a in abbrevs:
        d = os.path.join(inputs, a)
        s = json.load(open(os.path.join(d, "summary.json")))
        b = json.load(open(os.path.join(d, "exclusion_breakdown.json")))
        eb = b["excluded_breakdown"]
        L = lab.get(a, {})
        n_us = s["n_usable"]
        if n_us >= N_SUFFICIENT:
            suff = "yes"
        elif n_us >= N_MARGINAL:
            suff = "marginal"
        else:
            suff = "no"
        notes = []
        notes.append("%d unique haplotypes among %d usable records"
                     % (s["n_unique_usable"], n_us))
        if s["n_countries"] <= 1:
            notes.append("single country -> within-species pi confounded with "
                         "one local population")
        if s["host_available"] < n_us:
            notes.append("host missing for %d/%d" % (n_us - s["host_available"], n_us))
        if b["n_recoverable_full_orf"]:
            notes.append("%d further records carry a full-length N ORF but no CDS "
                         "feature (recoverable)" % b["n_recoverable_full_orf"])
        rows.append({
            "abbrev": a,
            "species": L.get("species", ""),
            "category": L.get("category", ""),
            "host_class": L.get("host_class", ""),
            "reference_accession": L.get("accession", ""),
            "esearch_records": s["esearch_count"],
            "records_M_or_L": eb.get("M_or_L_segment", 0),
            "records_partial_N": eb.get("N_cds_partial", 0),
            "records_no_cds_feature": (eb.get("S_no_cds_feature_partial", 0)
                                       + eb.get("S_no_cds_feature_FULL_ORF", 0)
                                       + eb.get("other_cds_only", 0)),
            "records_unclassifiable": eb.get("unclassifiable", 0),
            "excl_incomplete_cds": s["excl_partial_or_short"],
            "excl_ambiguity": s["excl_ambiguity"],
            "n_complete_n_cds": s["n_complete_n_cds"],
            "n_usable": n_us,
            "n_unique_usable": s["n_unique_usable"],
            "n_duplicate_records": s["n_duplicate_records"],
            "largest_duplicate_group": s["largest_duplicate_group"],
            "cds_len_min": s["cds_len_min"],
            "cds_len_modal": s["cds_len_modal"],
            "cds_len_max": s["cds_len_max"],
            "host_available": s["host_available"],
            "country_available": s["country_available"],
            "date_available": s["date_available"],
            "n_countries": s["n_countries"],
            "countries": ";".join(s["countries"]),
            "sufficient": suff,
            "notes": " | ".join(notes),
        })

    with open(out_tsv, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")

    # sanity: exclusion categories must account for every esearch hit
    for r in rows:
        tot = (r["records_M_or_L"] + r["records_partial_N"]
               + r["records_no_cds_feature"] + r["records_unclassifiable"]
               + r["excl_incomplete_cds"] + r["excl_ambiguity"] + r["n_usable"])
        flag = "OK" if tot == r["esearch_records"] else "MISMATCH"
        print("%-6s accounted %d / esearch %d  %s"
              % (r["abbrev"], tot, r["esearch_records"], flag))
    print("wrote", out_tsv)


if __name__ == "__main__":
    main()
