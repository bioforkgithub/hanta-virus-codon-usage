"""
CAI, ENC, mFE, and co-AI for the Kidney HKG (K-HKG) panel (Strauss et al. 2021,
18 candidate human kidney reference genes), computed identically to the Lung
HKG panel (02_codon_usage_metrics.py, 03_coai.py, 04_folding_energy.py), so
Andes virus ORFs can be benchmarked against both a pulmonary- and a
renal-tissue gene set (hantaviruses cause both HPS and HFRS/AKI).

Writes:
  results/kidney_hkg_metrics.tsv       (18 K-HKG genes, individually: CAI/ENC)
  results/mfe_kidney_hkg.tsv           (18 K-HKG genes, individually: mFE)
  results/codon_usage_kidney_hkg.tsv   (per-codon usage fraction, combined 18 K-HKG genes)
  results/coai_kidney.tsv              (co-AI for the K-HKG panel, appended alongside coai.tsv)
"""
import csv
import sys
from pathlib import Path
from scipy.stats import pearsonr

sys.path.insert(0, str(Path(__file__).parent))
from shared_utils import read_fasta, load_cai_weights, compute_cai, compute_enc, windowed_mfe, CODON_TABLE, SENSE_CODONS, strip_stop

BASE = Path(__file__).parent.parent
RESULTS = BASE / "results"
KIDNEY_FASTA = BASE.parent / "Kidney_HKG" / "kidney_hkg_cds.fasta"


def codon_usage_fraction(seqs):
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


def load_trna_counts():
    rows = list(csv.DictReader(open(RESULTS / "near_cognate_map.tsv"), delimiter="\t"))
    return {r["codon"]: int(r["final_tRNA_count_for_coAI"]) for r in rows}


def main():
    weights = load_cai_weights()
    kidney_recs = read_fasta(KIDNEY_FASTA)
    print(f"Loaded {len(kidney_recs)} kidney HKG CDS records")

    # --- CAI / ENC per gene ---
    rows = []
    for header, seq in kidney_recs:
        symbol = header.split("_")[0]
        cai = compute_cai(seq, weights)
        enc = compute_enc(seq)
        rows.append({"gene_symbol": symbol, "accession": header.split("_", 1)[1], "length_nt": len(strip_stop(seq)), "CAI": round(cai, 4), "ENC": round(enc, 3)})
    with open(RESULTS / "kidney_hkg_metrics.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} kidney HKG CAI/ENC rows")

    # --- mFE per gene ---
    mfe_rows = []
    for header, seq in kidney_recs:
        symbol = header.split("_")[0]
        seq_clean = strip_stop(seq)
        mfe = windowed_mfe(seq_clean)
        mfe_rows.append({"group": "Kidney_HKG", "name": symbol, "length_nt": len(seq_clean), "mFE_kcal_per_mol": round(mfe, 3)})
        print(f"Kidney HKG {symbol}: mFE = {mfe:.3f} kcal/mol")
    with open(RESULTS / "mfe_kidney_hkg.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mfe_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(mfe_rows)
    print(f"Wrote {len(mfe_rows)} kidney HKG mFE rows")

    # --- codon usage + co-AI ---
    write_codon_usage([seq for _, seq in kidney_recs], RESULTS / "codon_usage_kidney_hkg.tsv")
    trna_counts = load_trna_counts()
    usage_rows = list(csv.DictReader(open(RESULTS / "codon_usage_kidney_hkg.tsv"), delimiter="\t"))
    codons = [r["codon"] for r in usage_rows]
    usage = [float(r["usage_fraction"]) for r in usage_rows]
    trna = [trna_counts[c] for c in codons]
    r, p = pearsonr(usage, trna)
    with open(RESULTS / "coai_kidney.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["genome", "coAI_pearson_r", "p_value", "n_codons"])
        w.writerow(["Kidney_HKG_panel", round(r, 4), f"{p:.4g}", len(codons)])
    print(f"Kidney_HKG_panel: co-AI (Pearson r) = {r:.4f}, p = {p:.4g}, n = {len(codons)}")


if __name__ == "__main__":
    main()
