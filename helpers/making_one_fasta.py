import os
from Bio import SeqIO

# Define the root folder containing subfolders with fasta files
root_folder = "Hanta_All_data"
output_file = "All_hanta.fasta"

# List to store all sequences
all_sequences = []

# Traverse through all subfolders and fasta files
for subdir, dirs, files in os.walk(root_folder):
    for file in files:
        if file.endswith(".fna") or file.endswith(".fa"):
            file_path = os.path.join(subdir, file)
            # Parse fasta file and add sequences to the list
            with open(file_path, "r") as handle:
                sequences = list(SeqIO.parse(handle, "fasta"))
                all_sequences.extend(sequences)

# Write all sequences to the output fasta file
with open(output_file, "w") as output_handle:
    SeqIO.write(all_sequences, output_handle, "fasta")
