import os
"""
Download the 12 Human Lung Housekeeping Genes (L-HKG) used as the control
panel in the reference paper (mypaper.pdf, Table S2.1), by their exact
RefSeq accessions, and extract each gene's CDS.

Andes virus causes Hantavirus Pulmonary Syndrome (lung-tropic), so this
tissue-matched panel (decision D2) is directly reusable from the paper.
"""
import time
from pathlib import Path
from Bio import Entrez, SeqIO

Entrez.email = os.environ.get("ENTREZ_EMAIL", "bioforkgithub@outlook.com")

# Accession, gene symbol, as listed in mypaper.pdf Table S2.1
# NDUFA3 has no version in the paper; try latest current version via efetch by base accession.
ACCESSIONS = [
    ("NM_001101.5", "ACTB"),
    ("NM_001289746.2", "GAPDH"),
    ("NM_000291.4", "PGK1"),
    ("NM_012423.4", "RPL13A"),
    ("NM_053275.3", "RPLP0"),
    ("XM_005254549.3", "B2M"),
    ("NM_004168.4", "SDHA"),
    ("NM_000181.4", "GUSB"),
    ("NM_000190.4", "HMBS"),
    ("NM_000194.3", "HPRT1"),
    ("NM_003194.5", "TBP"),
    ("NM_004542", "NDUFA3"),
]

BASE = Path(__file__).parent
OUT_DIR = BASE.parent / "Lung_HKG"
OUT_DIR.mkdir(exist_ok=True)
(OUT_DIR / "gbk").mkdir(exist_ok=True)


def main():
    cds_fasta_lines = []
    report = []
    for acc, symbol in ACCESSIONS:
        h = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
        text = h.read()
        h.close()
        gbk_path = OUT_DIR / "gbk" / f"{symbol}_{acc}.gbk"
        gbk_path.write_text(text)

        rec = next(SeqIO.parse(gbk_path, "genbank"))
        cds_feats = [f for f in rec.features if f.type == "CDS"]
        if len(cds_feats) != 1:
            report.append((acc, symbol, rec.id, "UNEXPECTED_CDS_COUNT", len(cds_feats)))
            continue
        seq = str(cds_feats[0].extract(rec.seq)).upper()
        ok3 = len(seq) % 3 == 0
        ends_stop = seq[-3:] in {"TAA", "TAG", "TGA"}
        report.append((acc, symbol, rec.id, len(seq), f"x3={ok3} stop_end={ends_stop}"))
        cds_fasta_lines.append(f">{symbol}_{rec.id}\n{seq}\n")
        time.sleep(0.4)

    with open(OUT_DIR / "lung_hkg_cds.fasta", "w") as f:
        f.writelines(cds_fasta_lines)

    print(f"{'Requested':15s}{'Symbol':10s}{'Resolved':15s}{'Info'}")
    for acc, symbol, resolved, *rest in report:
        print(f"{acc:15s}{symbol:10s}{resolved:15s}{rest}")


if __name__ == "__main__":
    main()
