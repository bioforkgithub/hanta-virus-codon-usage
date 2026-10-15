#!/usr/bin/env python3
"""
Inventory complete nucleoprotein (S-segment) CDS per orthohantavirus species.

Phase 1 of the within-species piN/piS analysis. Downloads nothing into the
script directory: every species gets its own fresh directory under --outdir,
and all paths are taken as arguments. Run with `python3 -I` because the
downloaded data is untrusted.

Usage
-----
  python3 -I 01_inventory_N_cds.py \
      --outdir  .../PINPS_within_species_2026-10-12/inputs \
      --results .../PINPS_within_species_2026-10-12/results \
      --logs    .../PINPS_within_species_2026-10-12/logs \
      --group   environmental-A \
      --email   <ncbi contact email>

Determinism
-----------
No RNG is used in this phase. Records are sorted by accession before being
written, so the output files are byte-identical for a fixed GenBank state.
The GenBank query date is recorded in the log so the snapshot is identifiable.

Method
------
"Complete or near-complete N CDS" is defined on the PROTEIN database, not the
nucleotide database: a complete N CDS yields a full-length N protein, whereas a
partial CDS yields a short one. Orthohantavirus N is 428-436 aa (all six
species inventoried here translate to exactly 433 aa), and no other protein in
the orthohantavirus proteome falls in 410-450 aa (GPC ~1135 aa, L ~2150 aa), so
the length band alone is a sufficient and specific filter and no title keyword
match is needed. Title keyword matching was measured to MISS records (see
00_counts_via_eutils.md) and is therefore not used.

Each protein hit is mapped back to its nucleotide CDS via the /coded_by
qualifier, and that nucleotide span is what gets written to the per-species
FASTA. Curated duplicates (RefSeq, SwissProt) are flagged, not silently
dropped, because the task is to report duplication rather than hide it.
"""

import argparse
import hashlib
import json
import os
import re
import sys
import time
from collections import Counter

from Bio import Entrez, SeqIO

# --------------------------------------------------------------------------
# Species panel. Names are the ones used in the manuscript's labels.tsv; NCBI
# maps each to the current ICTV binomial itself (e.g. "Prospect Hill
# orthohantavirus" -> "Orthohantavirus prospectense"), which the log records.
# --------------------------------------------------------------------------
PANEL = {
    "environmental-A": [
        ("PHV",  "Prospect Hill orthohantavirus",  "U47136.1",    "environmental", "rodent"),
        ("TULV", "Tula orthohantavirus",           "NC_005227.1", "environmental", "rodent"),
        ("TPMV", "Thottapalayam orthohantavirus",  "NC_010704.1", "environmental", "shrew"),
        ("MJNV", "Imjin orthohantavirus",          "NC_034558.1", "environmental", "shrew"),
        ("NVAV", "Nova orthohantavirus",           "NC_034464.1", "environmental", "mole"),
        ("ASAV", "Asama orthohantavirus",          "NC_038273.1", "environmental", "mole"),
    ],
}

# Protein length band accepted as a complete / near-complete N.
AA_MIN, AA_MAX = 410, 450

# Ambiguity rules. The Epuyen work found N-runs masquerading as substitutions,
# so these are deliberately strict.
MAX_AMBIG_FRAC = 0.01   # >1% non-ACGT in the CDS -> exclude
MAX_AMBIG_RUN = 5       # any run of >5 consecutive Ns -> exclude

CODED_BY = re.compile(r"coded_by=\"?([A-Za-z0-9_.]+):(?:<)?(\d+)\.\.(?:>)?(\d+)")


def log(fh, msg):
    line = "%s\t%s" % (time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime()), msg)
    print(line)
    fh.write(line + "\n")
    fh.flush()


def esearch_count(term):
    h = Entrez.esearch(db="protein", term=term, retmax=0)
    rec = Entrez.read(h)
    h.close()
    return int(rec["Count"]), rec.get("QueryTranslation", "")


def esearch_ids(term, retmax=10000):
    h = Entrez.esearch(db="protein", term=term, retmax=retmax)
    rec = Entrez.read(h)
    h.close()
    return list(rec["IdList"]), rec.get("QueryTranslation", "")


def chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def ambiguity_profile(seq):
    """Return (n_ambiguous, ambiguous_fraction, longest_N_run)."""
    s = str(seq).upper()
    n_amb = sum(1 for c in s if c not in "ACGT")
    longest = 0
    run = 0
    for c in s:
        if c.upper() == "N":
            run += 1
            longest = max(longest, run)
        else:
            run = 0
    frac = (n_amb / len(s)) if s else 1.0
    return n_amb, frac, longest


def inventory_species(abbrev, sciname, ref_acc, category, host_class,
                      outdir, logfh):
    spdir = os.path.join(outdir, abbrev)
    os.makedirs(spdir, exist_ok=True)
    os.chmod(spdir, 0o700)

    org = '"%s"[Organism]' % sciname
    band = "%d:%d[SLEN]" % (AA_MIN, AA_MAX)

    # Tier 1: every protein record for the organism carrying an N title, any
    # length. Used only to report how many PARTIAL N records were excluded.
    t_titled, _ = esearch_count(
        "%s AND (nucleocapsid[Title] OR nucleoprotein[Title])" % org)
    # Tier 2: the complete/near-complete band (the operative definition).
    t_band, qtrans = esearch_count("%s AND %s" % (org, band))
    # Tier 3: band minus RefSeq copies.
    t_band_nonref, _ = esearch_count(
        "%s AND %s NOT srcdb_refseq[PROP]" % (org, band))
    # Tier 4: sanity check that nothing complete sits above the band.
    t_above, _ = esearch_count("%s AND %d:100000[SLEN] AND "
                               "(nucleocapsid[Title] OR nucleoprotein[Title])"
                               % (org, AA_MAX + 1))

    log(logfh, "%s\tquerytranslation\t%s" % (abbrev, qtrans))
    log(logfh, "%s\tcounts\ttitled_any_len=%d band=%d band_nonrefseq=%d "
               "above_band=%d" % (abbrev, t_titled, t_band, t_band_nonref,
                                  t_above))

    ids, _ = esearch_ids("%s AND %s" % (org, band))

    # --- protein summaries: length, source database, accession -------------
    prot = {}
    for batch in chunks(ids, 200):
        h = Entrez.esummary(db="protein", id=",".join(batch), retmode="xml")
        for d in Entrez.read(h):
            prot[str(d["AccessionVersion"])] = {
                "slen": int(d["Slen"]),
                "sourcedb": str(d.get("SourceDb", "")),
                "title": str(d.get("Title", "")),
            }
        h.close()
        time.sleep(0.4)

    # --- GenPept records: /coded_by -> nucleotide accession + CDS span -----
    for batch in chunks(sorted(prot), 100):
        h = Entrez.efetch(db="protein", id=",".join(batch),
                          rettype="gp", retmode="text")
        text = h.read()
        h.close()
        for blk in text.split("\n//\n"):
            m_acc = re.search(r"^VERSION\s+(\S+)", blk, re.M)
            m_cb = CODED_BY.search(blk)
            if not (m_acc and m_cb):
                continue
            acc = m_acc.group(1)
            if acc in prot:
                prot[acc]["nuc_acc"] = m_cb.group(1)
                prot[acc]["nuc_start"] = int(m_cb.group(2))
                prot[acc]["nuc_end"] = int(m_cb.group(3))
        time.sleep(0.4)

    # --- fetch the nucleotide records and slice out the CDS ----------------
    nuc_accs = sorted({v["nuc_acc"] for v in prot.values() if "nuc_acc" in v})
    raw_path = os.path.join(spdir, "nuccore_raw.gb")
    with open(raw_path, "w") as out:
        for batch in chunks(nuc_accs, 100):
            h = Entrez.efetch(db="nuccore", id=",".join(batch),
                              rettype="gb", retmode="text")
            out.write(h.read())
            h.close()
            time.sleep(0.4)
    os.chmod(raw_path, 0o600)

    meta = {}
    seqs = {}
    for rec in SeqIO.parse(raw_path, "genbank"):
        src = {}
        for f in rec.features:
            if f.type == "source":
                for q in ("host", "country", "geo_loc_name",
                          "isolation_source", "collection_date", "strain",
                          "isolate"):
                    if q in f.qualifiers:
                        src[q] = f.qualifiers[q][0]
        meta[rec.id.split(".")[0]] = src
        seqs[rec.id.split(".")[0]] = rec.seq

    # --- per-sequence QC table --------------------------------------------
    rows = []
    for pacc in sorted(prot):
        p = prot[pacc]
        nacc = p.get("nuc_acc")
        base = nacc.split(".")[0] if nacc else None
        if base is None or base not in seqs:
            rows.append({"protein_acc": pacc, "status": "no_coded_by",
                         "nuc_acc": nacc or ""})
            continue
        cds = seqs[base][p["nuc_start"] - 1:p["nuc_end"]]
        n_amb, frac, longest = ambiguity_profile(cds)
        src = meta.get(base, {})
        status = "usable"
        if p["slen"] < AA_MIN or p["slen"] > AA_MAX:
            status = "length_out_of_band"
        elif frac > MAX_AMBIG_FRAC:
            status = "excluded_ambiguous_frac"
        elif longest > MAX_AMBIG_RUN:
            status = "excluded_N_run"
        elif len(cds) % 3 != 0:
            status = "excluded_not_multiple_of_3"
        rows.append({
            "protein_acc": pacc,
            "nuc_acc": nacc,
            "sourcedb": p["sourcedb"],
            "is_curated_copy": p["sourcedb"] in ("refseq", "swiss_prot"),
            "aa_len": p["slen"],
            "cds_len": len(cds),
            "n_ambiguous": n_amb,
            "ambiguous_frac": round(frac, 6),
            "longest_N_run": longest,
            "md5": hashlib.md5(str(cds).upper().encode()).hexdigest(),
            "host": src.get("host", ""),
            "country": src.get("country", src.get("geo_loc_name", "")),
            "isolation_source": src.get("isolation_source", ""),
            "collection_date": src.get("collection_date", ""),
            "strain": src.get("strain", src.get("isolate", "")),
            "is_study_reference": (nacc == ref_acc),
            "status": status,
        })

    # --- duplicate accounting: report BOTH counts, never silently collapse -
    usable = [r for r in rows if r.get("status") == "usable"]
    primary = [r for r in usable if not r.get("is_curated_copy")]
    md5s = Counter(r["md5"] for r in primary)
    n_distinct = len(md5s)
    n_duplicate_records = sum(c - 1 for c in md5s.values() if c > 1)

    with open(os.path.join(spdir, "qc_table.tsv"), "w") as fh:
        cols = ["protein_acc", "nuc_acc", "sourcedb", "is_curated_copy",
                "aa_len", "cds_len", "n_ambiguous", "ambiguous_frac",
                "longest_N_run", "md5", "host", "country",
                "isolation_source", "collection_date", "strain",
                "is_study_reference", "status"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    os.chmod(os.path.join(spdir, "qc_table.tsv"), 0o600)

    # Population sample: keep duplicates (identical sequences are real data).
    fa = os.path.join(spdir, "%s_N_cds_all.fasta" % abbrev)
    with open(fa, "w") as fh:
        for r in sorted(primary, key=lambda x: x["nuc_acc"]):
            base = r["nuc_acc"].split(".")[0]
            p = prot[r["protein_acc"]]
            cds = seqs[base][p["nuc_start"] - 1:p["nuc_end"]]
            fh.write(">%s|%s\n%s\n" % (r["nuc_acc"], abbrev, str(cds).upper()))
    os.chmod(fa, 0o600)

    return {
        "abbrev": abbrev,
        "species": sciname,
        "category": category,
        "host_class": host_class,
        "reference_accession": ref_acc,
        "n_titled_any_length": t_titled,
        "n_complete_n_cds": t_band,
        "n_band_non_refseq": t_band_nonref,
        "n_above_band": t_above,
        "n_partial_excluded": t_titled - t_band,
        "n_curated_copies": sum(1 for r in rows if r.get("is_curated_copy")),
        "n_primary_records": len(primary),
        "n_usable": len(usable),
        "n_distinct_sequences": n_distinct,
        "n_duplicate_records": n_duplicate_records,
        "n_excluded_ambiguous": sum(
            1 for r in rows
            if str(r.get("status", "")).startswith("excluded")),
        "countries": len({r["country"] for r in primary if r["country"]}),
        "frac_with_host": (
            round(sum(1 for r in primary if r["host"]) / len(primary), 3)
            if primary else 0.0),
        "frac_with_country": (
            round(sum(1 for r in primary if r["country"]) / len(primary), 3)
            if primary else 0.0),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--logs", required=True)
    ap.add_argument("--group", default="environmental-A")
    ap.add_argument("--email", required=True)
    ap.add_argument("--api-key", default=None)
    a = ap.parse_args()

    Entrez.email = a.email
    if a.api_key:
        Entrez.api_key = a.api_key
    Entrez.tool = "PINPS_within_species_inventory"

    for d in (a.outdir, a.results, a.logs):
        os.makedirs(d, exist_ok=True)

    logpath = os.path.join(a.logs, "01_inventory_%s.log" % a.group)
    with open(logpath, "w") as logfh:
        log(logfh, "group=%s aa_band=%d:%d max_ambig_frac=%g max_N_run=%d"
                   % (a.group, AA_MIN, AA_MAX, MAX_AMBIG_FRAC, MAX_AMBIG_RUN))
        out = []
        for abbrev, sci, ref, cat, host in PANEL[a.group]:
            log(logfh, "=== %s (%s) ===" % (abbrev, sci))
            out.append(inventory_species(abbrev, sci, ref, cat, host,
                                         a.outdir, logfh))

        cols = list(out[0].keys())
        tsv = os.path.join(a.results, "inventory_%s.tsv" % a.group)
        with open(tsv, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in out:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        with open(os.path.join(a.results,
                               "inventory_%s.json" % a.group), "w") as fh:
            json.dump(out, fh, indent=2, sort_keys=True)
        log(logfh, "wrote %s" % tsv)


if __name__ == "__main__":
    sys.exit(main())
