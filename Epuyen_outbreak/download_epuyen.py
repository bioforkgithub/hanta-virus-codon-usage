import os
"""
Download the Epuyen 2018-2019 Andes virus person-to-person outbreak sequences.

28 patients, each with up to 3 genome segments (L=RdRp, M=glycoprotein
precursor Gn/Gc, S=nucleocapsid), all deposited to GenBank on 2020-09-30
(accessions MN258162-MN258256). Predecessor = Patient_1 (index case,
2018-11-03); successors = the remaining 27 patients.

Patient numbering follows the GenBank deposition, which runs 1-29 with no
Patient 21, giving 28 patients in total. Patient 2 has no L-segment
sequence deposited, and its M segment is partial (3345 nt, no stop codon),
so 80 rather than 81 patient-segment comparisons are evaluable downstream
(see pipeline/06_pnps_mutations.py, which skips absent segments).

Writes one GenBank flat file per accession to gbk/, and a manifest.tsv
mapping patient -> segment -> accession -> collection date.
"""
import re
import time
import csv
from pathlib import Path
from Bio import Entrez, SeqIO

Entrez.email = os.environ.get("ENTREZ_EMAIL", "bioforkgithub@outlook.com")

OUT_DIR = Path(__file__).parent
GBK_DIR = OUT_DIR / "gbk"
GBK_DIR.mkdir(exist_ok=True)

SEARCH_TERM = "Epuyen[All Fields] AND Andes[All Fields]"


def search_accessions():
    h = Entrez.esearch(db="nuccore", term=SEARCH_TERM, retmax=300)
    r = Entrez.read(h)
    return r["IdList"]


def download_all(ids):
    manifest = []
    batch = 20
    for i in range(0, len(ids), batch):
        chunk = ids[i : i + batch]
        h = Entrez.efetch(db="nuccore", id=",".join(chunk), rettype="gb", retmode="text")
        text = h.read()
        h.close()
        # split multi-record GenBank text on the "//" record terminator
        for rec_text in text.split("//\n"):
            rec_text = rec_text.strip()
            if not rec_text:
                continue
            rec_text += "\n//\n"
            m = re.search(r"^LOCUS\s+(\S+)", rec_text, re.M)
            if not m:
                continue
            accession = m.group(1)
            (GBK_DIR / f"{accession}.gbk").write_text(rec_text)
        time.sleep(0.4)

    for gbk_file in sorted(GBK_DIR.glob("*.gbk")):
        rec = next(SeqIO.parse(gbk_file, "genbank"))
        src = next(f for f in rec.features if f.type == "source")
        q = src.qualifiers
        isolate = q.get("isolate", ["NA"])[0]
        segment = q.get("segment", ["NA"])[0]
        pm = re.search(r"Patient_(\d+)_([\d/]+)", isolate)
        patient_num = pm.group(1) if pm else "NA"
        sample_date = pm.group(2) if pm else "NA"
        manifest.append(
            {
                "accession": rec.id,
                "isolate": isolate,
                "patient_num": patient_num,
                "sample_date": sample_date,
                "segment": segment,
                "length": len(rec.seq),
            }
        )

    manifest.sort(key=lambda r: (int(r["patient_num"]) if r["patient_num"].isdigit() else 999, r["segment"]))
    with open(OUT_DIR / "manifest.tsv", "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["accession", "isolate", "patient_num", "sample_date", "segment", "length"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(manifest)
    print(f"Wrote {len(manifest)} records to manifest.tsv")


if __name__ == "__main__":
    ids = search_accessions()
    print(f"Found {len(ids)} accessions")
    download_all(ids)
