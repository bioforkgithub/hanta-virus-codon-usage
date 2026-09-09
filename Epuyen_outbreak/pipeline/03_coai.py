"""
co-AI: Pearson correlation between codon usage (fraction) and human tRNA
gene copy number (near-cognate-mapped), across the 61 sense codons, per
genome - matching paper section 2.3/3.1 (Victor et al. tRNA co-adaptation
index).
"""
import csv
from pathlib import Path
from scipy.stats import pearsonr

BASE = Path(__file__).parent.parent
RESULTS = BASE / "results"


def load_trna_counts():
    rows = list(csv.DictReader(open(RESULTS / "near_cognate_map.tsv"), delimiter="\t"))
    return {r["codon"]: int(r["final_tRNA_count_for_coAI"]) for r in rows}


def coai_for_usage_table(usage_path, trna_counts):
    rows = list(csv.DictReader(open(usage_path), delimiter="\t"))
    codons = [r["codon"] for r in rows]
    usage = [float(r["usage_fraction"]) for r in rows]
    trna = [trna_counts[c] for c in codons]
    r, p = pearsonr(usage, trna)
    return r, p, len(codons)


def main():
    trna_counts = load_trna_counts()

    out_rows = []
    for label, path in [
        ("Andes_virus_predecessor_genome", RESULTS / "codon_usage_predecessor.tsv"),
        ("Lung_HKG_panel", RESULTS / "codon_usage_hkg.tsv"),
    ]:
        r, p, n = coai_for_usage_table(path, trna_counts)
        out_rows.append({"genome": label, "coAI_pearson_r": round(r, 4), "p_value": f"{p:.4g}", "n_codons": n})
        print(f"{label}: co-AI (Pearson r) = {r:.4f}, p = {p:.4g}, n = {n}")

    with open(RESULTS / "coai.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)


if __name__ == "__main__":
    main()
