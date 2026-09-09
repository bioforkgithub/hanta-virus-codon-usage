"""
Build a per-codon Homo sapiens cytosolic tRNA gene copy-number table from
the GtRNAdb hg38 high-confidence set already downloaded to hg38-tRNAs/
(hg38-tRNAs-confidence-set.out), matching the paper's methodology
(GtRNAdb 2.0, cognate tRNA copy number per anticodon -> codon).

Excludes:
  - SeC (selenocysteine) - not part of the standard 61-sense-codon table
  - 3 flagged anomalous entries (isotype mismatch x2, unexpected
    anticodon x1) - GtRNAdb itself flags these as unreliable
  - iMet (initiator tRNA) is counted separately from elongator Met, per
    the paper's convention (paper's Table S2.2 lists Met/CAT count=10,
    matching only the elongator "Met" isotype, not iMet)

Anticodon -> codon: codon = reverse complement of the anticodon (both
read 5'->3'), since codon:anticodon pairing is antiparallel.

Output: trna_codon_counts.tsv (codon, amino_acid, anticodon, tRNA_gene_count)
"""
import csv
from pathlib import Path
from Bio.Seq import Seq

BASE = Path(__file__).parent
SRC = BASE.parent / "hg38-tRNAs" / "hg38-tRNAs-confidence-set.out"
OUT = BASE / "trna_codon_counts.tsv"

BAD_NOTES = {"IPD:-95.80,isotype mismatch", "IPD:-43.10,isotype mismatch", "unexpected anticodon"}

# standard genetic code, DNA codon -> amino acid (3-letter, matching GtRNAdb isotype names)
CODON_TABLE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "*", "TAG": "*",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "*", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}


def main():
    lines = SRC.read_text().splitlines()[3:]  # skip 3 header lines
    counts = {}  # codon -> count
    anticodon_of_codon = {}
    skipped = {"SeC": 0, "iMet": 0, "bad_note": 0, "blank": 0}

    for line in lines:
        if not line.strip():
            skipped["blank"] += 1
            continue
        fields = line.split("\t")
        isotype = fields[4].strip()
        anticodon = fields[5].strip()
        note = fields[-1].strip()

        if isotype == "SeC":
            skipped["SeC"] += 1
            continue
        if isotype == "iMet":
            skipped["iMet"] += 1
            continue
        if note in BAD_NOTES:
            skipped["bad_note"] += 1
            continue

        codon = str(Seq(anticodon).reverse_complement())
        counts[codon] = counts.get(codon, 0) + 1
        anticodon_of_codon[codon] = anticodon

    rows = []
    for codon, aa in CODON_TABLE.items():
        if aa == "*":
            continue
        rows.append(
            {
                "codon": codon,
                "amino_acid": aa,
                "anticodon": anticodon_of_codon.get(codon, "NA"),
                "trna_gene_count": counts.get(codon, 0),
            }
        )
    rows.sort(key=lambda r: (r["amino_acid"], r["codon"]))

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["codon", "amino_acid", "anticodon", "trna_gene_count"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_zero = sum(1 for r in rows if r["trna_gene_count"] == 0)
    print(f"Wrote {len(rows)} sense codons to {OUT}")
    print(f"Codons with zero cognate tRNA (need near-cognate mapping): {n_zero}")
    print(f"Skipped: {skipped}")
    print(f"Total tRNA genes counted: {sum(counts.values())}")


if __name__ == "__main__":
    main()
