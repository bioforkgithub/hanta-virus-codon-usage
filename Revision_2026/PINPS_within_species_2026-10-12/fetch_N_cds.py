#!/usr/bin/env python3
import os  # noqa: E402  (added on publication, for ENTREZ_EMAIL)
"""
Inventory within-species complete nucleoprotein (S-segment) CDS in GenBank.

Phase 1 of the within-species piN/piS analysis. DOWNLOADS ONLY + CLASSIFIES.
Computes no piN/piS.

Usage:
    python3 -I fetch_N_cds.py <OUTROOT> <ABBREV> "<Organism name>" [--stage count|fetch]

OUTROOT is the untrusted-data root (inputs/). Each species gets its own fresh
subdirectory OUTROOT/<ABBREV>/.

Deterministic: no RNG used in this phase. Record order is NCBI-accession-sorted,
not download order, so reruns give byte-identical tables.
"""
import sys, os, time, json, re, urllib.parse, urllib.request, gzip, io

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "hantavirus_pinps_inventory"
EMAIL = os.environ.get("ENTREZ_EMAIL", "you@example.org")

# ---- search window -------------------------------------------------------
# Broad net: anything from a bare partial N fragment (300 nt) up to oversized
# records (5000 nt, e.g. S segment cloned in a vector / concatenated submission).
# Classification into complete vs partial happens downstream, so that the number
# of EXCLUDED partials is a measured quantity, not an artefact of the query.
SLEN_MIN, SLEN_MAX = 300, 5000

# ---- N CDS acceptance windows -------------------------------------------
# Reference N lengths (aa) from the labels.tsv reference accessions:
#   HTNV 429, SEOV 429, DOBV 429, PUUV 433, SNV 428  -> CDS 1287-1302 nt with stop
N_AA_MIN, N_AA_MAX = 380, 460          # candidate N CDS (excludes NSs ~70-90 aa)
N_NT_COMPLETE_MIN, N_NT_COMPLETE_MAX = 1260, 1320   # "complete" N CDS window
AMBIG_FRAC_MAX = 0.005                 # <=0.5% non-ACGT
MAX_N_RUN_MAX = 10                     # no N-run longer than 10 nt

N_PRODUCT_RE = re.compile(
    r"nucleocapsid|nucleoprotein|\bnucleo\b|\bN protein\b", re.I)
NSS_RE = re.compile(r"nonstructural|non-structural|\bNSs\b|NS protein", re.I)


def eget(endpoint, params, retries=5):
    params = dict(params)
    params.setdefault("tool", TOOL)
    params.setdefault("email", EMAIL)
    url = EUTILS + "/" + endpoint + "?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept-Encoding": "gzip"})
            with urllib.request.urlopen(req, timeout=180) as fh:
                raw = fh.read()
                if fh.headers.get("Content-Encoding") == "gzip":
                    raw = gzip.decompress(raw)
                return raw.decode("utf-8", "replace")
        except Exception as e:           # transient NCBI 429/500
            last = e
            time.sleep(3 * (attempt + 1))
    raise RuntimeError("eutils failed: %s (%s)" % (endpoint, last))


def build_query(organism):
    return ('("%s"[Organism]) AND %d:%d[SLEN] NOT patent[Properties]'
            % (organism, SLEN_MIN, SLEN_MAX))


def esearch_ids(organism):
    q = build_query(organism)
    txt = eget("esearch.fcgi", {"db": "nuccore", "term": q, "retmax": 0,
                               "usehistory": "n", "retmode": "json"})
    total = int(json.loads(txt)["esearchresult"]["count"])
    ids = []
    step = 5000
    for start in range(0, total, step):
        txt = eget("esearch.fcgi", {"db": "nuccore", "term": q,
                                    "retstart": start, "retmax": step,
                                    "retmode": "json", "sort": "accession"})
        ids.extend(json.loads(txt)["esearchresult"]["idlist"])
        time.sleep(0.4)
    return q, total, ids


def efetch_gb(ids, outpath, batch=150):
    """Download GenBank flatfiles into one concatenated file (untrusted data)."""
    with open(outpath, "w") as out:
        for i in range(0, len(ids), batch):
            chunk = ids[i:i + batch]
            txt = eget("efetch.fcgi", {"db": "nuccore", "id": ",".join(chunk),
                                       "rettype": "gbwithparts",
                                       "retmode": "text"})
            out.write(txt)
            if not txt.endswith("\n"):
                out.write("\n")
            time.sleep(0.4)
            sys.stderr.write("  fetched %d/%d\n" % (min(i + batch, len(ids)), len(ids)))


def max_run(seq, ch="N"):
    best = cur = 0
    for c in seq.upper():
        if c == ch:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def classify(gbfile, abbrev):
    from Bio import SeqIO
    from Bio.Seq import Seq
    rows = []
    for rec in SeqIO.parse(gbfile, "genbank"):
        src = {}
        for f in rec.features:
            if f.type == "source":
                src = f.qualifiers
                break
        meta = dict(
            host=(src.get("host") or [""])[0],
            country=(src.get("country") or src.get("geo_loc_name") or [""])[0],
            isolation_source=(src.get("isolation_source") or [""])[0],
            collection_date=(src.get("collection_date") or [""])[0],
            strain=(src.get("strain") or src.get("isolate") or [""])[0],
            organism=(src.get("organism") or [rec.annotations.get("organism", "")])[0],
        )
        defn = rec.description
        cds_found = []
        for f in rec.features:
            if f.type != "CDS":
                continue
            prod = " ".join(f.qualifiers.get("product", []))
            gene = " ".join(f.qualifiers.get("gene", []))
            note = " ".join(f.qualifiers.get("note", []))
            tag = prod + " | " + gene + " | " + note
            if NSS_RE.search(prod) or NSS_RE.search(gene):
                continue
            is_n = bool(N_PRODUCT_RE.search(tag)) or gene.strip().upper() in ("N", "S", "NP")
            try:
                nt = str(f.extract(rec.seq)).upper()
            except Exception:
                continue
            aa_len = len(nt) // 3
            if not is_n:
                # unannotated/ambiguously annotated CDS of exactly N size still counts
                if not (N_AA_MIN <= aa_len <= N_AA_MAX):
                    continue
            if not (N_AA_MIN <= aa_len <= N_AA_MAX):
                continue
            cds_found.append((f, nt, tag))

        if not cds_found:
            rows.append(dict(accession=rec.id, length=len(rec.seq), cds_nt=0,
                             source="no_N_CDS_annotated", definition=defn,
                             complete=False, usable=False, reason="no N CDS feature",
                             n_amb=0, amb_frac=0.0, max_N_run=0, cds_seq="", **meta))
            continue

        # keep the single best (longest) N candidate per record
        f, nt, tag = max(cds_found, key=lambda t: len(t[1]))
        fuzzy = ("<" in str(f.location)) or (">" in str(f.location))
        pseudo = "pseudo" in f.qualifiers or "pseudogene" in f.qualifiers
        partial_defn = bool(re.search(r"\bpartial\b", defn, re.I))
        n_amb = sum(1 for c in nt if c not in "ACGT")
        amb_frac = n_amb / len(nt) if nt else 1.0
        nrun = max_run(nt, "N")
        try:
            aa = str(Seq(nt).translate(table=1))
        except Exception:
            aa = ""
        starts_atg = nt.startswith("ATG")
        ends_stop = aa.endswith("*")
        internal_stop = "*" in aa[:-1]
        in_frame = len(nt) % 3 == 0

        reasons = []
        if fuzzy:
            reasons.append("fuzzy CDS boundary (<,>)")
        if pseudo:
            reasons.append("pseudo")
        if partial_defn and (fuzzy or not (N_NT_COMPLETE_MIN <= len(nt) <= N_NT_COMPLETE_MAX)):
            reasons.append("DEFINITION says partial")
        if not in_frame:
            reasons.append("CDS length not multiple of 3")
        if not (N_NT_COMPLETE_MIN <= len(nt) <= N_NT_COMPLETE_MAX):
            reasons.append("CDS length %d outside complete window %d-%d"
                           % (len(nt), N_NT_COMPLETE_MIN, N_NT_COMPLETE_MAX))
        if not starts_atg:
            reasons.append("no ATG start")
        if not ends_stop:
            reasons.append("no terminal stop")
        if internal_stop:
            reasons.append("internal stop codon")
        complete = not reasons

        qreasons = list(reasons)
        if amb_frac > AMBIG_FRAC_MAX:
            qreasons.append("ambiguity fraction %.4f > %.4f" % (amb_frac, AMBIG_FRAC_MAX))
        if nrun > MAX_N_RUN_MAX:
            qreasons.append("longest N-run %d > %d" % (nrun, MAX_N_RUN_MAX))
        usable = not qreasons

        rows.append(dict(accession=rec.id, length=len(rec.seq), cds_nt=len(nt),
                         source=tag.strip(), definition=defn,
                         complete=complete, usable=usable,
                         reason="; ".join(qreasons) if qreasons else "",
                         n_amb=n_amb, amb_frac=round(amb_frac, 6),
                         max_N_run=nrun, cds_seq=nt, **meta))
    rows.sort(key=lambda r: r["accession"])
    return rows


COLS = ["accession", "organism", "length", "cds_nt", "complete", "usable",
        "reason", "n_amb", "amb_frac", "max_N_run", "host", "country",
        "isolation_source", "collection_date", "strain", "source", "definition"]


def main():
    outroot, abbrev, organism = sys.argv[1], sys.argv[2], sys.argv[3]
    stage = "fetch"
    if "--stage" in sys.argv:
        stage = sys.argv[sys.argv.index("--stage") + 1]
    d = os.path.join(outroot, abbrev)
    os.makedirs(d, exist_ok=True)
    os.chmod(d, 0o700)

    q, total, ids = esearch_ids(organism)
    with open(os.path.join(d, "query.txt"), "w") as fh:
        fh.write(q + "\n")
        fh.write("esearch_count\t%d\n" % total)
        fh.write("uids_retrieved\t%d\n" % len(ids))
    sys.stderr.write("[%s] query hits=%d ids=%d\n" % (abbrev, total, len(ids)))
    if stage == "count":
        print(json.dumps({"abbrev": abbrev, "esearch_count": total}))
        return

    gb = os.path.join(d, "records.gb")
    if not (os.path.exists(gb) and os.path.getsize(gb) > 0):
        efetch_gb(ids, gb)

    rows = classify(gb, abbrev)

    with open(os.path.join(d, "records_table.tsv"), "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in COLS) + "\n")

    # FASTA of usable CDS (all of them, duplicates retained -- a population
    # sample's identical sequences are real data) + a de-duplicated FASTA.
    usable = [r for r in rows if r["usable"]]
    with open(os.path.join(d, "%s_N_cds_usable_all.fasta" % abbrev), "w") as fh:
        for r in usable:
            fh.write(">%s\n%s\n" % (r["accession"], r["cds_seq"]))
    seen = {}
    for r in usable:
        seen.setdefault(r["cds_seq"], []).append(r["accession"])
    with open(os.path.join(d, "%s_N_cds_usable_uniq.fasta" % abbrev), "w") as fh:
        for seq, accs in sorted(seen.items(), key=lambda kv: kv[1][0]):
            fh.write(">%s n_identical=%d\n%s\n" % (accs[0], len(accs), seq))
    with open(os.path.join(d, "duplicate_groups.tsv"), "w") as fh:
        fh.write("representative\tn_identical\tmembers\n")
        for seq, accs in sorted(seen.items(), key=lambda kv: kv[1][0]):
            fh.write("%s\t%d\t%s\n" % (accs[0], len(accs), ",".join(sorted(accs))))

    lens = sorted(r["cds_nt"] for r in usable)
    summ = dict(
        abbrev=abbrev, organism=organism, esearch_count=total,
        records_parsed=len(rows),
        no_N_cds=sum(1 for r in rows if r["cds_nt"] == 0),
        n_complete_n_cds=sum(1 for r in rows if r["complete"]),
        n_usable=len(usable),
        n_unique_usable=len(seen),
        n_duplicate_records=len(usable) - len(seen),
        largest_duplicate_group=(max(len(v) for v in seen.values()) if seen else 0),
        cds_len_min=(lens[0] if lens else 0), cds_len_max=(lens[-1] if lens else 0),
        cds_len_modal=(max(set(lens), key=lens.count) if lens else 0),
        host_available=sum(1 for r in usable if r["host"].strip()),
        country_available=sum(1 for r in usable if r["country"].strip()),
        isolsrc_available=sum(1 for r in usable if r["isolation_source"].strip()),
        date_available=sum(1 for r in usable if r["collection_date"].strip()),
        n_countries=len({r["country"].split(":")[0].strip() for r in usable
                         if r["country"].strip()}),
        countries=sorted({r["country"].split(":")[0].strip() for r in usable
                          if r["country"].strip()}),
        excl_partial_or_short=sum(1 for r in rows if r["cds_nt"] and not r["complete"]),
        excl_ambiguity=sum(1 for r in rows if r["complete"] and not r["usable"]),
    )
    with open(os.path.join(d, "summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2, sort_keys=True)
    print(json.dumps(summ, sort_keys=True))


if __name__ == "__main__":
    main()
