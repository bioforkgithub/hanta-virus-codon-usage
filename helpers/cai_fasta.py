import math
from Bio import SeqIO

# Codon usage weights (w) from the provided codon usage table (RSCU normalized to max per amino acid)
codon_weights = {
    "UUU": 0.46, "UUC": 0.54,
    "UCU": 0.19, "UCC": 0.22, "UCA": 0.15, "UCG": 0.05,
    "UAU": 0.44, "UAC": 0.56,
    "UGU": 0.46, "UGC": 0.54,
    "UGG": 1.00,
    "UUA": 0.08, "UUG": 0.13,
    "CUU": 0.13, "CUC": 0.20, "CUA": 0.07, "CUG": 0.40,
    "CCU": 0.29, "CCC": 0.32, "CCA": 0.28, "CCG": 0.11,
    "CAU": 0.42, "CAC": 0.58,
    "CAA": 0.27, "CAG": 0.73,
    "CGU": 0.08, "CGC": 0.18, "CGA": 0.11, "CGG": 0.20,
    "AUU": 0.36, "AUC": 0.47, "AUA": 0.17,
    "ACU": 0.25, "ACC": 0.36, "ACA": 0.28, "ACG": 0.11,
    "AAU": 0.47, "AAC": 0.53,
    "AAA": 0.43, "AAG": 0.57,
    "AGU": 0.15, "AGC": 0.24,
    "AGA": 0.21, "AGG": 0.21,
    "AUG": 1.00,
    "GUU": 0.18, "GUC": 0.24, "GUA": 0.12, "GUG": 0.46,
    "GCU": 0.27, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11,
    "GAU": 0.46, "GAC": 0.54,
    "GAA": 0.42, "GAG": 0.58,
    "GGU": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25
}

def dna_to_rna(dna_seq):
    """Convert a DNA sequence to RNA sequence by replacing T with U."""
    return dna_seq.upper().replace('T', 'U')

def calculate_cai(rna_sequence):
    """
    Calculate the Codon Adaptation Index (CAI) for a given RNA gene sequence.
    The sequence length must be a multiple of 3.
    """
    if len(rna_sequence) % 3 != 0:
        raise ValueError("Sequence length is not a multiple of 3.")

    weights = []
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if codon not in codon_weights:
            raise ValueError(f"Codon '{codon}' not found in codon weights table.")
        weights.append(codon_weights[codon])

    product = 1.0
    for w in weights:
        product *= w

    cai = product ** (1 / len(weights))
    return cai

def calculate_cai_from_fasta(fasta_file, output_file):
    """
    Read sequences from a FASTA file and calculate CAI for each.
    Writes results to the specified output file.
    Skips sequences with errors and logs the error messages.
    """
    with open(output_file, 'w') as out_f:
        out_f.write("Sequence_ID\tCAI_or_Error\n")
        for record in SeqIO.parse(fasta_file, "fasta"):
            dna_seq = str(record.seq)
            rna_seq = dna_to_rna(dna_seq)
            try:
                cai_value = calculate_cai(rna_seq)
                out_f.write(f"{record.id}\t{cai_value:.4f}\n")
            except ValueError as e:
                # Skip sequences that do not meet criteria, but log error
                out_f.write(f"{record.id}\tError: {str(e)}\n")

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python cai_fasta.py <input_fasta_file> <output_file>")
        sys.exit(1)

    fasta_path = sys.argv[1]
    output_path = sys.argv[2]

    calculate_cai_from_fasta(fasta_path, output_path)
    print(f"CAI calculation complete. Results saved to '{output_path}'.")
