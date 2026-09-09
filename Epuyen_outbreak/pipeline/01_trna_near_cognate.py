"""
Near-cognate tRNA mapping for codons with zero cognate tRNA genes.

Uses the SAME codon->anticodon near-cognate relationships established in
the reference paper's Table S2.2 (a general property of the human tRNA
pool / genetic code wobble rules - Homo sapiens has ~48 anticodons for 61
codons, so several codons rely on wobble decoding by a paralog tRNA gene
that primarily serves a different codon in the same amino acid family;
dos Reis et al. 2004, cited in the paper, is the standard reference for
these wobble/near-cognate rules). This mapping is independent of which
virus is being studied - it's a property of Homo sapiens, so it is
directly reusable here.

Applies it to OUR OWN tRNA gene counts (trna_codon_counts.tsv, built from
GtRNAdb hg38 independently in this project) rather than the paper's
published counts.
"""
import csv
from pathlib import Path

BASE = Path(__file__).parent.parent

# codon -> anticodon actually used for translation (own cognate anticodon
# if the codon has one; otherwise the near-cognate anticodon borrowed from
# a paralog tRNA within the same amino acid family), transcribed from the
# reference paper's Table S2.2.
NEAR_COGNATE_ANTICODON = {
    "GCA": "TGC", "GCC": "AGC", "GCG": "CGC", "GCT": "AGC",       # Ala
    "TGC": "GCA", "TGT": "GCA",                                    # Cys
    "GAC": "GTC", "GAT": "GTC",                                    # Asp
    "GAA": "TTC", "GAG": "CTC",                                    # Glu
    "TTC": "GAA", "TTT": "GAA",                                    # Phe
    "GGA": "TCC", "GGC": "GCC", "GGG": "CCC", "GGT": "GCC",       # Gly
    "CAC": "GTG", "CAT": "GTG",                                    # His
    "ATA": "TAT", "ATC": "GAT", "ATT": "AAT",                     # Ile
    "AAA": "TTT", "AAG": "CTT",                                    # Lys
    "CTA": "TAG", "CTC": "AAG", "CTG": "CAG", "CTT": "AAG",
    "TTA": "TAA", "TTG": "CAA",                                    # Leu
    "ATG": "CAT",                                                  # Met
    "AAC": "GTT", "AAT": "GTT",                                    # Asn
    "CCA": "TGG", "CCC": "AGG", "CCG": "CGG", "CCT": "AGG",       # Pro
    "CAA": "TTG", "CAG": "CTG",                                    # Gln
    "AGA": "TCT", "AGG": "CCT", "CGA": "TCG", "CGC": "ACG",
    "CGG": "CCG", "CGT": "ACG",                                    # Arg
    "AGC": "GCT", "AGT": "GCT", "TCA": "TGA", "TCC": "AGA",
    "TCG": "CGA", "TCT": "AGA",                                    # Ser
    "ACA": "TGT", "ACC": "AGT", "ACG": "CGT", "ACT": "AGT",       # Thr
    "GTA": "TAC", "GTC": "AAC", "GTG": "CAC", "GTT": "AAC",       # Val
    "TGG": "CCA",                                                  # Trp
    "TAC": "GTA", "TAT": "GTA",                                    # Tyr
}


def main():
    rows = list(csv.DictReader(open(BASE / "trna_codon_counts.tsv"), delimiter="\t"))
    anticodon_to_count = {r["anticodon"]: int(r["trna_gene_count"]) for r in rows if int(r["trna_gene_count"]) > 0}

    out_rows = []
    for r in rows:
        codon = r["codon"]
        own_count = int(r["trna_gene_count"])
        mapped_anticodon = NEAR_COGNATE_ANTICODON.get(codon)
        if mapped_anticodon is None:
            raise ValueError(f"No near-cognate mapping defined for codon {codon} - fix NEAR_COGNATE_ANTICODON")
        mapped_count = anticodon_to_count.get(mapped_anticodon)
        if mapped_count is None:
            raise ValueError(f"Mapped anticodon {mapped_anticodon} for codon {codon} has no gene count in our table")
        out_rows.append(
            {
                "codon": codon,
                "amino_acid": r["amino_acid"],
                "own_tRNA_count": own_count,
                "is_near_cognate": own_count == 0,
                "mapped_anticodon": mapped_anticodon,
                "final_tRNA_count_for_coAI": mapped_count,
            }
        )

    out_path = BASE / "results" / "near_cognate_map.tsv"
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    n_near = sum(1 for r in out_rows if r["is_near_cognate"])
    print(f"Wrote {len(out_rows)} codons to {out_path} ({n_near} near-cognate)")


if __name__ == "__main__":
    main()
