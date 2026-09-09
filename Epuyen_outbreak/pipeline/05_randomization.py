"""
Synonymous codon randomization (paper section 3.3 / S1.1): 100 randomized
sequences per predecessor ORF (sampling without replacement per amino
acid, preserving exact codon usage and protein sequence), mFE computed
for each, compared to the native ORF's mFE via a one-tailed paired t-test
(paper used t=-1.618, P<0.1, n=12 ORFs; we only have n=3 ORFs since
Hantavirus is segmented into far fewer genes than SARS-CoV-2 - this is a
real, expected difference in statistical power, not a bug. Documented as
a caveat.)
"""
import csv
import random
import sys
from pathlib import Path
from scipy.stats import ttest_rel

sys.path.insert(0, str(Path(__file__).parent))
from shared_utils import read_fasta, windowed_mfe, randomize_codons, strip_stop

BASE = Path(__file__).parent.parent
RESULTS = BASE / "results"
N_RANDOM = 100
SEED = 42


def main():
    rng = random.Random(SEED)
    per_orf_summary = []
    all_random_rows = []

    for seg in ["L", "M", "S"]:
        header, seq = read_fasta(BASE / "orfs" / "by_patient" / f"Patient_1_{seg}.fasta")[0]
        seq = strip_stop(seq)
        native_mfe = windowed_mfe(seq)

        random_mfes = []
        for i in range(N_RANDOM):
            rseq = randomize_codons(seq, rng)
            mfe = windowed_mfe(rseq)
            random_mfes.append(mfe)
            all_random_rows.append({"segment": seg, "replicate": i + 1, "mFE_kcal_per_mol": round(mfe, 3)})

        n_more_stable = sum(1 for m in random_mfes if m < native_mfe)  # more negative = more stable
        best_random = min(random_mfes)  # "stability maxima" among randomizations
        per_orf_summary.append(
            {
                "segment": seg,
                "native_mFE": round(native_mfe, 3),
                "randomized_mean_mFE": round(sum(random_mfes) / len(random_mfes), 3),
                "randomized_best_mFE": round(best_random, 3),
                "n_more_stable_than_native": n_more_stable,
                "pct_more_stable_than_native": round(100 * n_more_stable / N_RANDOM, 1),
            }
        )
        print(f"{seg}: native={native_mfe:.3f}, {n_more_stable}/{N_RANDOM} randomizations more stable, best={best_random:.3f}")

    with open(RESULTS / "randomization_summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(per_orf_summary[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(per_orf_summary)

    with open(RESULTS / "randomization_all_replicates.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(all_random_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(all_random_rows)

    natives = [r["native_mFE"] for r in per_orf_summary]
    best_randoms = [r["randomized_best_mFE"] for r in per_orf_summary]
    t_stat, p_two_tailed = ttest_rel(natives, best_randoms)
    p_one_tailed = p_two_tailed / 2  # one-tailed, matching paper's reporting convention
    print(f"\nPaired one-tailed t-test (native vs. best-of-100-randomized mFE, n=3 ORFs): t={t_stat:.3f}, p(one-tailed)={p_one_tailed:.4f}")
    print("CAVEAT: n=3 (Hantavirus has 3 ORFs vs SARS-CoV-2's 12) - this test is severely underpowered; interpret with caution.")

    with open(RESULTS / "randomization_ttest.tsv", "w") as f:
        f.write("statistic\tvalue\n")
        f.write(f"t\t{t_stat:.4f}\n")
        f.write(f"p_one_tailed\t{p_one_tailed:.4f}\n")
        f.write(f"n_orfs\t{len(natives)}\n")
        f.write("caveat\tn=3 ORFs (vs paper's n=12) - severely underpowered, interpret with caution\n")


if __name__ == "__main__":
    main()
