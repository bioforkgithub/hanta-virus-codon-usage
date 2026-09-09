from Bio import SeqIO
import csv

def read_ensembl_ids_from_csv(file_path, id_column, delimiter=','):
    ensembl_ids = set()
    with open(file_path, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=delimiter)
        for row in reader:
            ensembl_id = row[id_column].strip()
            if ensembl_id:
                ensembl_ids.add(ensembl_id)
    return ensembl_ids

def extract_fasta_sequences_biopython(fasta_file, include_ids, exclude_ids, output_file):
    with open(output_file, 'w') as f_out:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            ensembl_id = record.id.split('.')[0]
            if ensembl_id in include_ids and ensembl_id not in exclude_ids:
                SeqIO.write(record, f_out, 'fasta')

# Paths to input files
fasta_file = 'Homo_sapiens.GRCh38.cds.all.fa'
housekeeping_csv = 'Housekeeping_GenesHuman.csv'
moststable_csv = 'MostStableExpressed.csv'

# Read Ensembl IDs from Housekeeping_GenesHuman.csv (comma-separated)
hk_ensembl_ids = read_ensembl_ids_from_csv(housekeeping_csv, id_column='Ensembl', delimiter=';')

# Read Ensembl IDs from MostStableExpressed.csv (tab-separated)
stable_ensembl_ids = read_ensembl_ids_from_csv(moststable_csv, id_column='Ensembl ID', delimiter='\t')

# Extract sequences for housekeeping genes excluding those in most stable expressed genes
extract_fasta_sequences_biopython(fasta_file, hk_ensembl_ids, stable_ensembl_ids, 'hk_human.fasta')

# Extract sequences for most stable expressed genes
extract_fasta_sequences_biopython(fasta_file, stable_ensembl_ids, set(), 'stableHK.fasta')
