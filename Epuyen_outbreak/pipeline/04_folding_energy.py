"""
Windowed minimum folding energy (mFE) for Patient_1 (predecessor)'s 3 ORFs
and the 12 Lung HKG mRNAs, matching paper section 2.4/S1.2.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from shared_utils import read_fasta, windowed_mfe, strip_stop

BASE = Path(__file__).parent.parent
RESULTS = BASE / "results"


def main():
    rows = []

    for seg in ["L", "M", "S"]:
        header, seq = read_fasta(BASE / "orfs" / "by_patient" / f"Patient_1_{seg}.fasta")[0]
        seq = strip_stop(seq)
        mfe = windowed_mfe(seq)
        product = " ".join(header.split()[2:])
        rows.append({"group": "Andes_virus_predecessor", "name": f"{seg}_{product}", "length_nt": len(seq), "mFE_kcal_per_mol": round(mfe, 3)})
        print(f"Predecessor {seg} ({product}): mFE = {mfe:.3f} kcal/mol")

    hkg_recs = read_fasta(BASE.parent / "Lung_HKG" / "lung_hkg_cds.fasta")
    for header, seq in hkg_recs:
        symbol = header.split("_")[0]
        seq = strip_stop(seq)
        mfe = windowed_mfe(seq)
        rows.append({"group": "Lung_HKG", "name": symbol, "length_nt": len(seq), "mFE_kcal_per_mol": round(mfe, 3)})
        print(f"HKG {symbol}: mFE = {mfe:.3f} kcal/mol")

    with open(RESULTS / "mfe_predecessor_hkg.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} mFE rows to results/mfe_predecessor_hkg.tsv")


if __name__ == "__main__":
    main()
