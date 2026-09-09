"""
Extract clean, in-frame ORF nucleotide sequences from the Epuyen outbreak
GenBank records (gbk/*.gbk), using the single annotated CDS feature per
segment (L=RdRp, M=glycoprotein precursor, S=nucleoprotein).

Outputs (all under orfs/):
  by_segment/{L,M,S}_all_patients.fasta  - one multi-fasta per segment,
      one record per patient (for EMBOSS cusp/chips per-ORF-type analyses,
      and for per-ORF CAI/folding-energy tables analogous to paper Table S2.3)
  by_patient/Patient{N}_{segment}.fasta  - one file per patient per segment
  by_patient_combined/Patient{N}_allORFs.fasta - all available ORFs for one
      patient concatenated as separate FASTA records (for whole-genome-level
      EMBOSS cusp/chips Nc, analogous to paper Table 1's per-strain ENC)

Also writes orfs/orf_summary.tsv with length/stop-codon sanity checks.
"""
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

BASE = Path(__file__).parent
GBK_DIR = BASE / "gbk"
ORF_DIR = BASE / "orfs"
(ORF_DIR / "by_segment").mkdir(parents=True, exist_ok=True)
(ORF_DIR / "by_patient").mkdir(parents=True, exist_ok=True)
(ORF_DIR / "by_patient_combined").mkdir(parents=True, exist_ok=True)

STOP_CODONS = {"TAA", "TAG", "TGA"}


def load_manifest():
    with open(BASE / "manifest.tsv") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def extract_cds(gbk_path):
    rec = next(SeqIO.parse(gbk_path, "genbank"))
    cds_feats = [f for f in rec.features if f.type == "CDS"]
    if len(cds_feats) != 1:
        raise ValueError(f"{gbk_path}: expected exactly 1 CDS, found {len(cds_feats)}")
    feat = cds_feats[0]
    seq = feat.extract(rec.seq)
    product = feat.qualifiers.get("product", ["NA"])[0]
    return str(seq).upper(), product


def main():
    manifest = load_manifest()
    by_segment_records = {"L": [], "M": [], "S": []}
    by_patient = {}
    summary_rows = []

    for row in manifest:
        acc = row["accession"]
        seg = row["segment"]
        pnum = row["patient_num"]
        if seg not in ("L", "M", "S"):
            continue
        gbk_path = GBK_DIR / f"{acc.split('.')[0]}.gbk"
        seq, product = extract_cds(gbk_path)

        length_ok = len(seq) % 3 == 0
        internal_stop = any(
            seq[i : i + 3] in STOP_CODONS for i in range(0, len(seq) - 3, 3)
        )
        ends_in_stop = seq[-3:] in STOP_CODONS

        summary_rows.append(
            {
                "patient_num": pnum,
                "segment": seg,
                "accession": acc,
                "product": product,
                "length": len(seq),
                "multiple_of_3": length_ok,
                "internal_stop_before_end": internal_stop,
                "ends_in_stop": ends_in_stop,
            }
        )

        label = f"Patient_{pnum}"
        header = f"{label}_{seg} {acc} {product}"

        # per-patient, per-segment file
        pf = ORF_DIR / "by_patient" / f"{label}_{seg}.fasta"
        pf.write_text(f">{header}\n{seq}\n")

        by_segment_records[seg].append((header, seq))
        by_patient.setdefault(pnum, {})[seg] = (header, seq)

    # by_segment multi-fasta
    for seg, records in by_segment_records.items():
        out = ORF_DIR / "by_segment" / f"{seg}_all_patients.fasta"
        with open(out, "w") as f:
            for header, seq in records:
                f.write(f">{header}\n{seq}\n")

    # by_patient_combined (all available ORFs for that patient)
    for pnum, segs in by_patient.items():
        out = ORF_DIR / "by_patient_combined" / f"Patient_{pnum}_allORFs.fasta"
        with open(out, "w") as f:
            for seg in ("L", "M", "S"):
                if seg in segs:
                    header, seq = segs[seg]
                    f.write(f">{header}\n{seq}\n")

    with open(ORF_DIR / "orf_summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "patient_num",
                "segment",
                "accession",
                "product",
                "length",
                "multiple_of_3",
                "internal_stop_before_end",
                "ends_in_stop",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(summary_rows)

    n_bad = sum(
        1 for r in summary_rows if not r["multiple_of_3"] or r["internal_stop_before_end"]
    )
    print(f"Extracted {len(summary_rows)} ORFs across {len(by_patient)} patients.")
    print(f"Sequences failing sanity checks (not x3 or internal stop): {n_bad}")
    print(f"Sequences ending in a stop codon: {sum(1 for r in summary_rows if r['ends_in_stop'])}")


if __name__ == "__main__":
    main()
