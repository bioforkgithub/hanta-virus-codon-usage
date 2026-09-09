import os
import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = os.environ.get("ENTREZ_EMAIL", "bioforkgithub@outlook.com")

tsv_file = "hanta_virus.tsv"
gbk_folder = "gbk_file"
genomes_folder = "genomes"
metadata_file = "metadata.tsv"

# Ensure directories exist before any file operations
os.makedirs(gbk_folder, exist_ok=True)
os.makedirs(genomes_folder, exist_ok=True)

df = pd.read_csv(tsv_file, sep='\t')

def download_gbk(accession, out_folder):
    try:
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
        data = handle.read()
        handle.close()
        if not data.strip():
            print(f"Warning: Empty GBK data for {accession}")
            return None
        gbk_path = os.path.join(out_folder, f"{accession}.gbk")
        with open(gbk_path, "w") as f:
            f.write(data)
        return gbk_path
    except Exception as e:
        print(f"Error downloading {accession}: {e}")
        return None

def extract_orfs(gbk_path, out_folder):
    try:
        records = list(SeqIO.parse(gbk_path, "genbank"))
        if not records:
            print(f"No records found in GBK file {gbk_path} for ORF extraction")
            return
        for idx, record in enumerate(records):
            orfs = []
            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    orf_seq = feature.qualifiers["translation"][0]
                    orf_id = feature.qualifiers.get("protein_id", [f"unknown_{len(orfs)+1}"])[0]
                    orf_desc = feature.qualifiers.get("product", ["unknown_product"])[0]
                    orfs.append((orf_id, orf_desc, orf_seq))
            if not orfs:
                print(f"No CDS features with translation found in record {idx+1} of {gbk_path}")
                continue
            fasta_name = os.path.splitext(os.path.basename(gbk_path))[0]
            if len(records) > 1:
                fasta_name += f"_record{idx+1}"
            fasta_path = os.path.join(out_folder, f"{fasta_name}_orfs.fasta")
            with open(fasta_path, "w") as fasta_file:
                for orf_id, orf_desc, orf_seq in orfs:
                    fasta_file.write(f">{orf_id} {orf_desc}\n{orf_seq}\n")
    except Exception as e:
        print(f"Error extracting ORFs from {gbk_path}: {e}")

metadata_rows = []
for accession in df['Assembly Accession']:
    gbk_path = download_gbk(accession, gbk_folder)
    if gbk_path:
        extract_orfs(gbk_path, genomes_folder)
        try:
            records = list(SeqIO.parse(gbk_path, "genbank"))
            if not records:
                print(f"No records found for metadata extraction in {gbk_path}")
                continue
            # Aggregate metadata from all records or use first record as representative
            # Here using first record for simplicity
            first_record = records[0]
            metadata_rows.append({
                "Assembly Accession": accession,
                "Organism": first_record.annotations.get("organism", ""),
                "Sequence_Length": len(first_record.seq),
                "Description": first_record.description
            })
        except Exception as e:
            print(f"Error reading metadata from {gbk_path}: {e}")

# Ensure gbk_folder still exists before saving metadata (in case of external deletion)
if not os.path.exists(gbk_folder):
    os.makedirs(gbk_folder, exist_ok=True)

metadata_df = pd.DataFrame(metadata_rows)
metadata_df.to_csv(os.path.join(gbk_folder, metadata_file), sep='\t', index=False)
