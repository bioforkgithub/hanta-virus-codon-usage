#!/usr/bin/env python3
"""Prepare the 22-species S-segment codon alignment for ML tree + PAML.

Analyses covered here: R2.2 (tree-based dN/dS), R3.2
(model-based estimate, group ordering), R4.3 (phylogenetic codon models),
R4.12 (ML tree replacing NJ/p-distance).

Inputs are read from Work1/Genus_Diversity and never modified.
Outputs: clean-named FASTA, PAML sequential PHYLIP, and a label table.
"""
import sys, csv
from pathlib import Path
from Bio import SeqIO

GD = Path(sys.argv[1])       # Work1/Genus_Diversity
OUT = Path(sys.argv[2])      # analysis inputs dir

aln_in = GD / "genus_diversity_codon_aligned.fasta"
man_in = GD / "genus_diversity_manifest.tsv"

# --- labels -----------------------------------------------------------------
HOST_OF = {  # host class parsed from the manifest 'note' field, checked by hand
    "rodent": {"HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV",
               "CHOV", "LANV", "PHV", "TULV"},
    "shrew":  {"TPMV", "MJNV", "RKPV", "SWSV", "CBNV", "JJUV"},
    "mole":   {"NVAV", "ASAV", "OXBV"},
    "bat":    {"QZNV"},
}
host_lookup = {ab: h for h, s in HOST_OF.items() for ab in s}

rows = list(csv.DictReader(man_in.open(), delimiter="\t"))
labels = {}
for r in rows:
    ab = r["abbrev"]
    labels[ab] = {
        "abbrev": ab,
        "species": r["species"],
        "accession": r["accession"],
        "category": r["category"],
        "host_class": host_lookup[ab],
        "is_rodent": str(host_lookup[ab] == "rodent"),
        "note": r["note"],
    }

# --- alignment --------------------------------------------------------------
recs = list(SeqIO.parse(aln_in, "fasta"))
clean = []
for rec in recs:
    ab = rec.id.split("|")[0]
    if ab not in labels:
        sys.exit(f"alignment name {rec.id!r} not in manifest")
    rec.id = ab
    rec.description = ""
    clean.append(rec)

n = len(clean)
L = len(clean[0].seq)
if any(len(r.seq) != L for r in clean):
    sys.exit("alignment is not flush")
if L % 3:
    sys.exit(f"alignment length {L} is not a multiple of 3")

# codon-level sanity: internal stop codons break codeml
STOPS = {"TAA", "TAG", "TGA"}
problems = []
for rec in clean:
    s = str(rec.seq).upper()
    cods = [s[i:i + 3] for i in range(0, L, 3)]
    internal = [(j + 1, c) for j, c in enumerate(cods[:-1]) if c in STOPS]
    if internal:
        problems.append((rec.id, internal))

OUT.mkdir(parents=True, exist_ok=True)
SeqIO.write(clean, OUT / "S22_codon_aln.fasta", "fasta")

# PAML sequential phylip
with (OUT / "S22_codon_aln.phy").open("w") as fh:
    fh.write(f" {n} {L}\n")
    for rec in clean:
        fh.write(f"{rec.id}\n{str(rec.seq).upper()}\n")

with (OUT / "labels.tsv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
        "abbrev", "species", "accession", "category", "host_class",
        "is_rodent", "note"])
    w.writeheader()
    for ab in [r.id for r in clean]:
        w.writerow(labels[ab])

# --- report -----------------------------------------------------------------
npath = sum(1 for v in labels.values() if v["category"] == "pathogenic")
nenv = sum(1 for v in labels.values() if v["category"] == "environmental")
print(f"sequences        : {n}")
print(f"alignment length : {L} nt = {L // 3} codons")
print(f"pathogenic       : {npath}")
print(f"environmental    : {nenv}")
print("host classes     : " + ", ".join(
    f"{h}={sum(1 for v in labels.values() if v['host_class'] == h)}"
    for h in ("rodent", "shrew", "mole", "bat")))
print("rodent-borne environmental (R3.1 confounding check): " + ", ".join(
    sorted(ab for ab, v in labels.items()
           if v["category"] == "environmental" and v["host_class"] == "rodent")))
if problems:
    print("\nWARNING internal stop codons (would break codeml):")
    for ab, lst in problems:
        print(f"  {ab}: {lst[:5]}{' …' if len(lst) > 5 else ''}")
else:
    print("internal stop codons: none")
