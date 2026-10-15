#!/usr/bin/env python3
"""Phase 2: inventory complete nucleoprotein (S-segment) CDS per species.

Run with:
  python3 -I 02_inventory_environmental_B.py <inputs_dir> <results_dir> <logs_dir>

Reads the UNTRUSTED GenBank flatfiles under <inputs_dir>/<ABBREV>/records.gb
(parsed as data only -- nothing from that tree is imported or executed).

Deterministic: no RNG; every list is sorted by accession. Seed is irrelevant in
this phase and none is used (the resampling phase will state its own seed).
"""
import json
import os
import re
import sys
from collections import Counter

from Bio import SeqIO

ABBREVS = ["OXBV", "RKPV", "SWSV", "CBNV", "QZNV", "JJUV"]

# --- N CDS recognition -------------------------------------------------------
N_PRODUCT = re.compile(
    r"nucleocapsid|nucleoprotein|\bprotein\s+N\b|\bN\s+protein\b", re.I)
NOT_N_PRODUCT = re.compile(
    r"nonstructural|non-structural|\bNSs?\b|glycoprotein|polymerase|"
    r"\bGPC\b|\bGn\b|\bGc\b|\bL\s+protein\b|RNA-dependent", re.I)

# --- completeness / quality thresholds (stated in the notes file) ------------
LEN_MIN_COMPLETE = 1260   # nt, incl. stop; orthohantavirus N is 428-433 aa
LEN_MAX_COMPLETE = 1320
LEN_MIN_NEAR = 1200       # 'near-complete' floor
MAX_NONACGT_FRAC = 0.01   # >1% ambiguous -> drop
MAX_N_RUN = 10            # Epuyen work: N-runs masquerade as substitutions

# --- sufficiency rule --------------------------------------------------------
MIN_USABLE = 10
MIN_HAPLOTYPES = 4


def longest_n_run(seq):
    runs = re.findall(r"N+", seq.upper())
    return max((len(r) for r in runs), default=0)


def is_partial(feature):
    return "<" in str(feature.location) or ">" in str(feature.location)


def find_orf(seq, min_len=LEN_MIN_COMPLETE, max_len=LEN_MAX_COMPLETE):
    """Longest ATG..stop ORF on the + strand within the length window."""
    s = str(seq).upper()
    stops = {"TAA", "TAG", "TGA"}
    best = None
    for frame in range(3):
        i = frame
        while i < len(s) - 2:
            if s[i:i + 3] == "ATG":
                j = i
                while j < len(s) - 2:
                    if s[j:j + 3] in stops:
                        L = j + 3 - i
                        if min_len <= L <= max_len and (best is None or L > best[1] - best[0]):
                            best = (i, j + 3)
                        break
                    j += 3
            i += 3
    return best


def src_meta(rec):
    for f in rec.features:
        if f.type == "source":
            q = f.qualifiers
            return {
                "host": (q.get("host") or [""])[0],
                "country": (q.get("country") or q.get("geo_loc_name") or [""])[0],
                "isolation_source": (q.get("isolation_source") or [""])[0],
                "collection_date": (q.get("collection_date") or [""])[0],
                "organism": (q.get("organism") or [""])[0],
                "segment": (q.get("segment") or [""])[0],
            }
    return {k: "" for k in ("host", "country", "isolation_source",
                            "collection_date", "organism", "segment")}


def classify_record(rec):
    """Return (list_of_candidate_dicts, record_level_reason_if_none)."""
    meta = src_meta(rec)
    defn = rec.description or ""
    out = []
    cds_seen = 0
    for f in rec.features:
        if f.type != "CDS":
            continue
        cds_seen += 1
        prod = " ".join(f.qualifiers.get("product", []) +
                        f.qualifiers.get("gene", []) +
                        f.qualifiers.get("note", []))
        gene = " ".join(f.qualifiers.get("gene", []))
        is_n = bool(N_PRODUCT.search(prod)) or gene.strip().upper() in ("N", "S")
        if NOT_N_PRODUCT.search(prod) and not N_PRODUCT.search(prod):
            is_n = False
        if not is_n:
            continue
        try:
            nt = str(f.extract(rec.seq)).upper()
        except Exception:
            continue
        out.append({"source": "cds_feature", "nt": nt,
                    "partial_location": is_partial(f),
                    "codon_start": (f.qualifiers.get("codon_start") or ["1"])[0],
                    "product": prod[:120]})
    if out:
        return out, meta, defn, cds_seen, None

    # no annotated N CDS: try ORF rescue if this looks like an S-segment record
    looks_s = (meta["segment"].strip().upper() in ("S", "SEGMENT S") or
               re.search(r"\bS\s+segment\b|nucleocapsid|nucleoprotein", defn, re.I))
    if looks_s and len(rec.seq) >= LEN_MIN_COMPLETE:
        hit = find_orf(rec.seq)
        if hit:
            nt = str(rec.seq[hit[0]:hit[1]]).upper()
            return ([{"source": "orf_rescue", "nt": nt,
                      "partial_location": False, "codon_start": "1",
                      "product": "longest-ORF (no CDS annotation)"}],
                    meta, defn, cds_seen, None)
        return [], meta, defn, cds_seen, "S-like, no CDS and no in-window ORF"
    return [], meta, defn, cds_seen, (
        "no N CDS annotated and record too short / not S-segment")


def main():
    inputs_dir, results_dir, logs_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    per_seq_rows, summary, log = [], {}, []

    for ab in ABBREVS:
        gb = os.path.join(inputs_dir, ab, "records.gb")
        recs = list(SeqIO.parse(gb, "genbank"))
        recs.sort(key=lambda r: r.id)

        n_records = len(recs)
        excl = Counter()
        cand_rows = []

        for rec in recs:
            cands, meta, defn, cds_seen, reason = classify_record(rec)
            if not cands:
                excl["no_N_CDS_in_record"] += 1
                per_seq_rows.append({
                    "abbrev": ab, "accession": rec.id, "rec_len": len(rec.seq),
                    "cds_source": "", "cds_len": "", "status": "excluded",
                    "exclusion_reason": reason or "no N CDS",
                    "n_nonacgt": "", "max_n_run": "",
                    "host": meta["host"], "country": meta["country"],
                    "isolation_source": meta["isolation_source"],
                    "collection_date": meta["collection_date"],
                    "definition": defn[:140]})
                continue
            # if several N CDS in one record, keep the longest
            cand = max(cands, key=lambda c: len(c["nt"]))
            nt = cand["nt"]
            L = len(nt)
            nonacgt = sum(1 for c in nt if c not in "ACGT")
            nrun = longest_n_run(nt)

            status, reason = "usable", ""
            if cand["partial_location"] or cand["codon_start"] != "1":
                status, reason = "excluded", "partial CDS (fuzzy location / codon_start>1)"
            elif L < LEN_MIN_NEAR:
                status, reason = "excluded", f"too short ({L} nt < {LEN_MIN_NEAR})"
            elif L < LEN_MIN_COMPLETE:
                status, reason = "excluded", f"near-complete only ({L} nt)"
            elif L > LEN_MAX_COMPLETE:
                status, reason = "excluded", f"over-long ({L} nt > {LEN_MAX_COMPLETE})"
            elif L % 3 != 0:
                status, reason = "excluded", f"length not multiple of 3 ({L})"
            elif not nt.startswith("ATG"):
                status, reason = "excluded", "no ATG start"
            elif nonacgt / L > MAX_NONACGT_FRAC:
                status, reason = "excluded", f"ambiguous {nonacgt}/{L} > 1%"
            elif nrun > MAX_N_RUN:
                status, reason = "excluded", f"N-run {nrun} > {MAX_N_RUN}"
            if status == "excluded":
                excl[reason.split(" (")[0]] += 1

            row = {"abbrev": ab, "accession": rec.id, "rec_len": len(rec.seq),
                   "cds_source": cand["source"], "cds_len": L,
                   "status": status, "exclusion_reason": reason,
                   "n_nonacgt": nonacgt, "max_n_run": nrun,
                   "host": meta["host"], "country": meta["country"],
                   "isolation_source": meta["isolation_source"],
                   "collection_date": meta["collection_date"],
                   "definition": defn[:140]}
            per_seq_rows.append(row)
            if status == "usable":
                cand_rows.append((rec.id, nt, meta, cand["source"]))

        # ---- complete-CDS pool (before quality filter) ----
        complete_pool = [r for r in per_seq_rows
                         if r["abbrev"] == ab and r["cds_len"] != ""
                         and LEN_MIN_COMPLETE <= r["cds_len"] <= LEN_MAX_COMPLETE
                         and not str(r["exclusion_reason"]).startswith("partial")]
        usable = sorted(cand_rows, key=lambda t: t[0])
        n_usable = len(usable)

        # ---- de-duplication (recorded, not silently collapsed) ----
        hap = Counter(nt for _, nt, _, _ in usable)
        n_hap = len(hap)
        n_dupes = n_usable - n_hap
        hap_sizes = sorted(hap.values(), reverse=True)

        # ---- metadata availability ----
        host_vals = [m["host"] for _, _, m, _ in usable if m["host"]]
        ctry_vals = [m["country"] for _, _, m, _ in usable if m["country"]]
        isos_vals = [m["isolation_source"] for _, _, m, _ in usable
                     if m["isolation_source"]]
        countries = sorted({c.split(":")[0].strip() for c in ctry_vals})
        hosts = sorted(set(host_vals))

        lens = sorted(r["cds_len"] for r in per_seq_rows
                      if r["abbrev"] == ab and r["cds_len"] != "")
        orf_rescued = sum(1 for _, _, _, s in usable if s == "orf_rescue")

        sufficient = (n_usable >= MIN_USABLE and n_hap >= MIN_HAPLOTYPES)

        # write the per-species FASTA next to its own untrusted data
        fa = os.path.join(inputs_dir, ab, "N_cds_usable.fasta")
        with open(fa, "w") as fh:
            for acc, nt, m, s in usable:
                fh.write(f">{acc}|{ab}|len{len(nt)}|host={m['host'] or 'NA'}"
                         f"|country={m['country'] or 'NA'}|{s}\n{nt}\n")

        summary[ab] = {
            "n_genbank_records": n_records,
            "n_complete_N_cds": len(complete_pool),
            "n_usable": n_usable,
            "n_distinct_haplotypes": n_hap,
            "n_duplicate_sequences": n_dupes,
            "haplotype_size_spectrum": hap_sizes,
            "n_orf_rescued": orf_rescued,
            "cds_len_min": lens[0] if lens else None,
            "cds_len_max": lens[-1] if lens else None,
            "cds_len_median": lens[len(lens) // 2] if lens else None,
            "exclusions": dict(sorted(excl.items())),
            "frac_with_host": round(len(host_vals) / n_usable, 3) if n_usable else None,
            "frac_with_country": round(len(ctry_vals) / n_usable, 3) if n_usable else None,
            "frac_with_isolation_source": round(len(isos_vals) / n_usable, 3) if n_usable else None,
            "n_countries": len(countries), "countries": countries,
            "hosts": hosts,
            "sufficient": sufficient,
            "fasta": fa,
        }
        log.append(f"{ab}: records={n_records} complete={len(complete_pool)} "
                   f"usable={n_usable} haplotypes={n_hap} dupes={n_dupes} "
                   f"orf_rescued={orf_rescued} sufficient={sufficient}")
        print(log[-1], flush=True)

    # ---- outputs ----
    tsv = os.path.join(results_dir, "inventory_environmental-B.tsv")
    cols = ["abbrev", "category", "host_class", "n_genbank_records",
            "n_complete_N_cds", "n_usable", "n_distinct_haplotypes",
            "n_duplicate_sequences", "n_orf_rescued", "cds_len_min",
            "cds_len_median", "cds_len_max", "frac_with_host",
            "frac_with_country", "frac_with_isolation_source", "n_countries",
            "countries", "hosts", "sufficient"]
    HOSTCLASS = {"OXBV": "mole", "RKPV": "shrew", "SWSV": "shrew",
                 "CBNV": "shrew", "QZNV": "bat", "JJUV": "shrew"}
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for ab in ABBREVS:
            s = summary[ab]
            fh.write("\t".join(str(x) for x in [
                ab, "environmental", HOSTCLASS[ab], s["n_genbank_records"],
                s["n_complete_N_cds"], s["n_usable"], s["n_distinct_haplotypes"],
                s["n_duplicate_sequences"], s["n_orf_rescued"], s["cds_len_min"],
                s["cds_len_median"], s["cds_len_max"], s["frac_with_host"],
                s["frac_with_country"], s["frac_with_isolation_source"],
                s["n_countries"], ";".join(s["countries"]) or "NA",
                ";".join(s["hosts"]) or "NA", s["sufficient"]]) + "\n")

    det = os.path.join(results_dir, "per_sequence_environmental-B.tsv")
    dcols = ["abbrev", "accession", "rec_len", "cds_source", "cds_len", "status",
             "exclusion_reason", "n_nonacgt", "max_n_run", "host", "country",
             "isolation_source", "collection_date", "definition"]
    with open(det, "w") as fh:
        fh.write("\t".join(dcols) + "\n")
        for r in sorted(per_seq_rows, key=lambda r: (r["abbrev"], r["accession"])):
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in dcols) + "\n")

    with open(os.path.join(results_dir, "inventory_environmental-B.json"), "w") as fh:
        json.dump({"thresholds": {
            "LEN_MIN_COMPLETE": LEN_MIN_COMPLETE, "LEN_MAX_COMPLETE": LEN_MAX_COMPLETE,
            "LEN_MIN_NEAR": LEN_MIN_NEAR, "MAX_NONACGT_FRAC": MAX_NONACGT_FRAC,
            "MAX_N_RUN": MAX_N_RUN, "MIN_USABLE": MIN_USABLE,
            "MIN_HAPLOTYPES": MIN_HAPLOTYPES},
            "species": summary}, fh, indent=2, sort_keys=True)

    with open(os.path.join(logs_dir, "02_inventory.log"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("wrote", tsv)


if __name__ == "__main__":
    main()
