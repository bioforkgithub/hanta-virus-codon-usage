"""
CAI, ENC, and codon-usage tables for:
  - each of the 28 patients' combined genome (available ORFs concatenated)
  - each patient's individual ORFs (L/M/S) - per-ORF CAI/ENC/mFE table analogous to paper Table S2.3
  - the 12 Lung HKG genes

Writes:
  results/strain_genome_metrics.tsv   (per patient, combined-genome ENC/CAI)
  results/predecessor_orf_metrics.tsv (Patient_1's 3 ORFs, individually)
  results/hkg_metrics.tsv             (12 HKG genes, individually)
  results/codon_usage_predecessor.tsv (per-codon usage fraction, Patient_1 combined genome - for co-AI)
  results/codon_usage_hkg.tsv         (per-codon usage fraction, combined 12 HKG genes - for co-AI)
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from shared_utils import read_fasta, load_cai_weights, compute_cai, compute_enc, CODON_TABLE, SENSE_CODONS, strip_stop

BASE = Path(__file__).parent.parent
RESULTS = BASE / "results"
RESULTS.mkdir(exist_ok=True)


def codon_usage_fraction(seqs):
    """seqs: list of nucleotide strings. Returns {codon: fraction of all sense codons used}."""
    counts = {c: 0 for c in SENSE_CODONS}
    total = 0
    for seq in seqs:
        seq = strip_stop(seq)
        for i in range(0, len(seq) - len(seq) % 3, 3):
            c = seq[i : i + 3]
            if c in counts:
                counts[c] += 1
                total += 1
    return {c: (counts[c] / total if total else 0.0) for c in SENSE_CODONS}


def write_codon_usage(seqs, out_path):
    frac = codon_usage_fraction(seqs)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["codon", "amino_acid", "usage_fraction"])
        for c in sorted(frac):
            w.writerow([c, CODON_TABLE[c], f"{frac[c]:.6f}"])


def main():
    weights = load_cai_weights()

    manifest = list(csv.DictReader(open(BASE / "manifest.tsv"), delimiter="\t"))
    patients = sorted(set(r["patient_num"] for r in manifest), key=int)

    # --- per-patient combined-genome metrics ---
    strain_rows = []
    predecessor_seqs = None
    for p in patients:
        combined_path = BASE / "orfs" / "by_patient_combined" / f"Patient_{p}_allORFs.fasta"
        if not combined_path.exists():
            continue
        recs = read_fasta(combined_path)
        seqs = [seq for _, seq in recs]
        segs_present = "".join(sorted(set(h.split()[0].rsplit("_", 1)[1] for h, _ in recs)))
        combined_seq = "".join(strip_stop(s) for s in seqs)
        try:
            cai = compute_cai(combined_seq, weights)
        except ValueError as e:
            cai = None
        enc = compute_enc(combined_seq)
        strain_rows.append(
            {
                "patient_num": p,
                "role": "predecessor" if p == "1" else "successor",
                "segments_present": segs_present,
                "n_segments": len(recs),
                "combined_length_nt": len(combined_seq),
                "ENC": round(enc, 3),
                "CAI": round(cai, 4) if cai else "NA",
            }
        )
        if p == "1":
            predecessor_seqs = recs  # [(header, seq), ...] for L, M, S

    with open(RESULTS / "strain_genome_metrics.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(strain_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(strain_rows)
    print(f"Wrote {len(strain_rows)} strain genome metrics")

    # --- predecessor per-ORF metrics (Table S2.3 analog) ---
    orf_rows = []
    for header, seq in predecessor_seqs:
        seg = header.split()[0].rsplit("_", 1)[1]
        product = " ".join(header.split()[2:])
        cai = compute_cai(seq, weights)
        enc = compute_enc(seq)
        orf_rows.append({"segment": seg, "product": product, "length_nt": len(strip_stop(seq)), "CAI": round(cai, 4), "ENC": round(enc, 3)})
    with open(RESULTS / "predecessor_orf_metrics.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(orf_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(orf_rows)
    print(f"Wrote {len(orf_rows)} predecessor ORF metrics")

    # --- HKG metrics ---
    hkg_recs = read_fasta(BASE.parent / "Lung_HKG" / "lung_hkg_cds.fasta")
    hkg_rows = []
    for header, seq in hkg_recs:
        symbol = header.split("_")[0]
        cai = compute_cai(seq, weights)
        enc = compute_enc(seq)
        hkg_rows.append({"gene_symbol": symbol, "accession": header.split("_", 1)[1], "length_nt": len(strip_stop(seq)), "CAI": round(cai, 4), "ENC": round(enc, 3)})
    with open(RESULTS / "hkg_metrics.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(hkg_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(hkg_rows)
    print(f"Wrote {len(hkg_rows)} HKG metrics")

    # --- codon usage tables for co-AI ---
    write_codon_usage([seq for _, seq in predecessor_seqs], RESULTS / "codon_usage_predecessor.tsv")
    write_codon_usage([seq for _, seq in hkg_recs], RESULTS / "codon_usage_hkg.tsv")
    print("Wrote codon usage tables for co-AI")


if __name__ == "__main__":
    main()
