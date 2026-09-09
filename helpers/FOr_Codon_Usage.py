####################################################################################################################
def check_and_write_multiples_of_three(fasta_file):
    from Bio import SeqIO
    from pathlib import Path
    input_path = Path(fasta_file)
    output_file = input_path.with_name(input_path.stem + "multiples_of_three" + input_path.suffix)

    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if len(record.seq) % 3 == 0:
                SeqIO.write(record, out_f, "fasta")
####################################################################################################################