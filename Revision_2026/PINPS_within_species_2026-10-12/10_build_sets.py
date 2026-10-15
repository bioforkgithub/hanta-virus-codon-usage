#!/usr/bin/env python3
"""
Phase 2 step 1: assemble the per-species nucleoprotein CDS sequence sets for
within-species piN/piS.

Reads ONLY from inputs/<ABBREV>/ (untrusted downloaded GenBank data) and writes
to results/sets/.  Computes no diversity statistics.

Two input flavours, both already on disk from phase 1:
  flavour "table" -- inputs/<A>/records_table.tsv written by scripts/fetch_N_cds.py
                     (columns accession, organism, usable, cds_seq is NOT in the
                     table, so the CDS is re-read from <A>_N_cds_usable_all.fasta)
  flavour "fasta" -- inputs/<A>/N_cds_refined.fasta written by the
                     environmental-B phase-1 scripts, header carries org=/host=/country=

Sets emitted per species (all are row-subsets of the same usable pool, so one
alignment serves all of them):
  lato  -- every usable complete N CDS under the [Organism] subtree
  ss    -- sensu stricto: lato minus records whose ORGANISM string names a
           different, separately-named virus lineage nested under that species
  (uniq / thin / animal subsets are defined downstream in 12_pinps.py from the
   metadata table written here, not as separate FASTAs)

Usage:  python3 -I 10_build_sets.py <ANALYSIS_ROOT>
Deterministic: no RNG.  Output order is accession-sorted.
"""
import sys, os, csv, re, json

# abbrev, category, host_class, flavour, regex of FOREIGN lineage names to drop
# for sensu stricto.  Every token below was taken from the observed ORGANISM
# strings in inputs/<A>/records_table.tsv (printed in logs/10_build_sets.log).
SPECIES = [
    ("HTNV", "pathogenic",   "rodent", "table", r"\bAmur\b|\bSoochong"),
    ("SEOV", "pathogenic",   "rodent", "table", r"\bGou|Sapporo rat"),
    ("DOBV", "pathogenic",   "rodent", "table", r"\bKurkino\b"),
    ("PUUV", "pathogenic",   "rodent", "table", r"\bHokkaido\b|\bMuju\b|Puumala-like"),
    ("SNV",  "pathogenic",   "rodent", "table", r"\bNew York\b"),
    ("ANDV", "pathogenic",   "rodent", "table",
     r"\bTunari\b|\bBermejo\b|\bLechiguanas\b|\bOran\b|\bMaciel\b|\bPergamino\b"
     r"|\bNeembucu\b|Araraquara"),
    ("CHOV", "pathogenic",   "rodent", "table", None),
    ("TULV", "environmental", "rodent", "table", r"\bAdler\b"),
    ("NVAV", "environmental", "mole",   "table", None),
    ("MJNV", "environmental", "shrew",  "table", None),
    ("TPMV", "environmental", "shrew",  "table", None),
    ("PHV",  "environmental", "rodent", "table", None),
    ("ASAV", "environmental", "mole",   "table", None),
    ("SWSV", "environmental", "shrew",  "fasta", None),
    ("CBNV", "environmental", "shrew",  "fasta", r"\bLianghe\b"),
    ("JJUV", "environmental", "shrew",  "fasta", None),
    ("RKPV", "environmental", "shrew",  "fasta", None),
    ("OXBV", "environmental", "mole",   "fasta", None),
    ("QZNV", "environmental", "bat",    "fasta", None),
]

HUMAN_RE = re.compile(r"homo sapiens|\bhuman\b|patient", re.I)


def read_fasta(path):
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def load_table(indir, abbrev):
    tsv = os.path.join(indir, "records_table.tsv")
    fa = os.path.join(indir, "%s_N_cds_usable_all.fasta" % abbrev)
    seqs = {h.split()[0]: s.upper() for h, s in read_fasta(fa)}
    recs = []
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["usable"] != "True":
                continue
            acc = r["accession"]
            if acc not in seqs:
                raise SystemExit("usable record %s missing from %s" % (acc, fa))
            recs.append(dict(
                accession=acc, organism=r["organism"], seq=seqs[acc],
                host=r["host"], country=r["country"],
                collection_date=r["collection_date"], strain=r["strain"]))
    return recs


def load_refined_fasta(indir):
    fa = os.path.join(indir, "N_cds_refined.fasta")
    recs = []
    for h, s in read_fasta(fa):
        parts = h.split("|")
        kv = {}
        for p in parts[3:]:
            if "=" in p:
                k, v = p.split("=", 1)
                kv[k] = v
        recs.append(dict(
            accession=parts[0], organism=kv.get("org", ""), seq=s.upper(),
            host=("" if kv.get("host", "NA") == "NA" else kv.get("host", "")),
            country=("" if kv.get("country", "NA") == "NA" else kv.get("country", "")),
            collection_date="", strain=(parts[2] if len(parts) > 2 else "")))
    return recs


def main():
    root = sys.argv[1]
    indir_root = os.path.join(root, "inputs")
    setdir = os.path.join(root, "results", "sets")
    os.makedirs(setdir, exist_ok=True)
    os.chmod(setdir, 0o700)

    meta_rows = []
    summary = {}
    for abbrev, cat, host_class, flavour, foreign in SPECIES:
        indir = os.path.join(indir_root, abbrev)
        recs = load_table(indir, abbrev) if flavour == "table" else load_refined_fasta(indir)
        recs.sort(key=lambda r: r["accession"])
        fre = re.compile(foreign) if foreign else None
        orgs, dropped = {}, []
        for r in recs:
            r["is_foreign"] = bool(fre and fre.search(r["organism"]))
            r["is_human"] = bool(HUMAN_RE.search(r["host"]))
            r["country_top"] = r["country"].split(":")[0].strip()
            r["year"] = (re.search(r"(\d{4})", r["collection_date"]).group(1)
                         if re.search(r"(\d{4})", r["collection_date"]) else "")
            orgs[r["organism"]] = orgs.get(r["organism"], 0) + 1
            if r["is_foreign"]:
                dropped.append(r["accession"] + " [" + r["organism"] + "]")
        ss = [r for r in recs if not r["is_foreign"]]

        # FASTA for alignment: the full lato pool (ss is a row subset of it)
        with open(os.path.join(setdir, "%s.lato.fasta" % abbrev), "w") as fh:
            for r in recs:
                fh.write(">%s\n%s\n" % (r["accession"], r["seq"]))

        for r in recs:
            meta_rows.append(dict(abbrev=abbrev, category=cat, host_class=host_class,
                                  accession=r["accession"], organism=r["organism"],
                                  in_sensu_stricto=(not r["is_foreign"]),
                                  cds_nt=len(r["seq"]), host=r["host"],
                                  is_human=r["is_human"], country=r["country_top"],
                                  year=r["year"], strain=r["strain"]))
        summary[abbrev] = dict(
            category=cat, host_class=host_class, flavour=flavour,
            n_lato=len(recs), n_sensu_stricto=len(ss),
            n_dropped_foreign=len(recs) - len(ss),
            dropped_foreign=dropped,
            organism_labels={k: v for k, v in sorted(orgs.items(), key=lambda kv: -kv[1])},
            n_human_ss=sum(1 for r in ss if r["is_human"]),
            n_nonhuman_ss=sum(1 for r in ss if r["host"].strip() and not r["is_human"]),
            n_host_missing_ss=sum(1 for r in ss if not r["host"].strip()),
            countries_ss=sorted({r["country_top"] for r in ss if r["country_top"]}),
            n_country_year_cells_ss=len({(r["country_top"], r["year"]) for r in ss
                                         if r["country_top"] and r["year"]}),
            cds_len_set=sorted({len(r["seq"]) for r in recs}),
        )
        print("%-5s lato=%-4d ss=%-4d dropped=%-3d human=%-4d nonhuman=%-4d countries=%d"
              % (abbrev, len(recs), len(ss), len(recs) - len(ss),
                 summary[abbrev]["n_human_ss"], summary[abbrev]["n_nonhuman_ss"],
                 len(summary[abbrev]["countries_ss"])))

    cols = ["abbrev", "category", "host_class", "accession", "organism",
            "in_sensu_stricto", "cds_nt", "host", "is_human", "country", "year", "strain"]
    with open(os.path.join(root, "results", "sets", "sequence_metadata.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in meta_rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")
    with open(os.path.join(root, "results", "sets", "set_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
