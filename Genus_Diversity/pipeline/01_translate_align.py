"""
Protein-guided codon alignment for the 22-species genus-wide nucleoprotein
(S-segment) dataset. Because these species are far more divergent than the
near-identical Andes virus outbreak genomes, a direct nucleotide alignment
is not reliable for codon-based dN/dS - instead:

  1. Translate each CDS to protein.
  2. Align the 22 protein sequences with MAFFT (L-INS-i via --auto).
  3. Back-translate the protein alignment to a codon-level nucleotide
     alignment (each aligned residue column maps to its original codon;
     gap columns get "---").

Output: genus_diversity_protein_aligned.fasta, genus_diversity_codon_aligned.fasta
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "Epuyen_outbreak" / "pipeline"))
from shared_utils import read_fasta, CODON_TABLE, strip_stop

BASE = Path(__file__).parent.parent
import os
MAFFT = os.environ.get("MAFFT", "mafft")

# shared_utils.CODON_TABLE stores 3-letter amino acid codes (e.g. "Leu"); MAFFT
# needs single-letter protein sequences.
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}


def translate(seq):
    seq = strip_stop(seq)
    codons = [seq[i : i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    return "".join(AA3_TO_1.get(CODON_TABLE.get(c), "X") for c in codons)


def main():
    recs = read_fasta(BASE / "genus_diversity_S_cds.fasta")
    protein_path = BASE / "genus_diversity_protein.fasta"
    with open(protein_path, "w") as f:
        for header, seq in recs:
            prot = translate(seq)
            f.write(f">{header}\n{prot}\n")
    print(f"Wrote {len(recs)} translated protein sequences")

    aligned_path = BASE / "genus_diversity_protein_aligned.fasta"
    with open(aligned_path, "w") as out:
        result = subprocess.run([MAFFT, "--auto", "--quiet", str(protein_path)], stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print("MAFFT STDERR:", result.stderr)
        sys.exit(1)
    print(f"Wrote protein alignment to {aligned_path.name}")

    # back-translate: map each aligned protein sequence's residues back to
    # their original codons, preserving gap columns.
    nt_lookup = {header: strip_stop(seq) for header, seq in recs}
    aligned_prot = read_fasta(aligned_path)

    codon_aligned_lines = []
    for header, aln_prot in aligned_prot:
        nt_seq = nt_lookup[header]
        codons = [nt_seq[i : i + 3] for i in range(0, len(nt_seq) - len(nt_seq) % 3, 3)]
        out_codons = []
        ci = 0
        for aa in aln_prot:
            if aa == "-":
                out_codons.append("---")
            else:
                out_codons.append(codons[ci])
                ci += 1
        assert ci == len(codons), f"{header}: consumed {ci} codons, expected {len(codons)}"
        codon_aligned_lines.append(f">{header}\n{''.join(out_codons)}\n")

    codon_aligned_path = BASE / "genus_diversity_codon_aligned.fasta"
    with open(codon_aligned_path, "w") as f:
        f.writelines(codon_aligned_lines)
    print(f"Wrote codon-level alignment to {codon_aligned_path.name}")
    print(f"Alignment length: {len(aligned_prot[0][1])} aa columns / {len(aligned_prot[0][1]) * 3} nt columns")


if __name__ == "__main__":
    main()
