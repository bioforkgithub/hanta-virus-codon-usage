"""
pn/ps polymorphism + mutation classification for each of the 26 successor
patients vs. Patient_1 (predecessor), per ORF (L/M/S) - matching paper
section 2.4/3.4/3.6, Tables 2-4.

Implements the standard Nei & Gojobori (1986) counting method:
  - synonymous/nonsynonymous SITES per codon (fraction of the 3 possible
    single-nt substitutions at each position that are synonymous)
  - synonymous/nonsynonymous DIFFERENCES between two aligned sequences,
    using pathway-averaging for codons differing at >1 position
  pn = Nd/N (nonsynonymous differences per nonsynonymous site)
  ps = Sd/S (synonymous differences per synonymous site)

Also classifies each individual codon change: synonymous/nonsynonymous,
amino acid polarity change, and tRNA gene-copy-number impact of the
change (cognate/near-cognate categories, paper Table 2 style).
"""
import csv
import sys
from itertools import permutations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from shared_utils import read_fasta, strip_stop, CODON_TABLE, POLARITY

BASE = Path(__file__).parent.parent
RESULTS = BASE / "results"
BASES = "ACGT"


def synonymous_sites_per_codon(codon):
    """Return S_i (fractional synonymous sites) for one codon, 0-3."""
    aa = CODON_TABLE[codon]
    s = 0.0
    for pos in range(3):
        syn_count = 0
        for b in BASES:
            if b == codon[pos]:
                continue
            mutant = codon[:pos] + b + codon[pos + 1 :]
            if CODON_TABLE.get(mutant) == aa:
                syn_count += 1
        s += syn_count / 3
    return s


def gene_sites(seq):
    """Total S (synonymous) and N (nonsynonymous) sites over a CDS."""
    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    S = sum(synonymous_sites_per_codon(c) for c in codons if CODON_TABLE.get(c, "*") != "*")
    N = 3 * len([c for c in codons if CODON_TABLE.get(c, "*") != "*"]) - S
    return S, N


def codon_pathway_diff(codon1, codon2):
    """
    Nei-Gojobori pathway-averaged synonymous/nonsynonymous difference
    count between two codons (0, 1, 2, or 3 nt differences). Returns
    (Sd_contribution, Nd_contribution) as floats (averaged over valid
    pathways that don't pass through a stop codon; if all pathways pass
    through a stop, falls back to unweighted average over all pathways).
    """
    diff_positions = [i for i in range(3) if codon1[i] != codon2[i]]
    if not diff_positions:
        return 0.0, 0.0
    if len(diff_positions) == 1:
        pos = diff_positions[0]
        aa1, aa2 = CODON_TABLE.get(codon1), CODON_TABLE.get(codon2)
        return (1.0, 0.0) if aa1 == aa2 else (0.0, 1.0)

    # multiple differences: enumerate all orderings of the differing positions
    pathways = []
    for order in permutations(diff_positions):
        current = codon1
        syn, nonsyn = 0, 0
        valid = True
        for pos in order:
            nxt = current[:pos] + codon2[pos] + current[pos + 1 :]
            aa_cur = CODON_TABLE.get(current)
            aa_nxt = CODON_TABLE.get(nxt)
            if aa_nxt == "*":
                valid = False
            if aa_cur == aa_nxt:
                syn += 1
            else:
                nonsyn += 1
            current = nxt
        pathways.append((syn, nonsyn, valid))

    valid_pathways = [(s, n) for s, n, v in pathways if v] or [(s, n) for s, n, v in pathways]
    Sd = sum(s for s, n in valid_pathways) / len(valid_pathways)
    Nd = sum(n for s, n in valid_pathways) / len(valid_pathways)
    return Sd, Nd


def pn_ps(ref_seq, alt_seq):
    length = min(len(ref_seq), len(alt_seq))
    length -= length % 3
    ref_codons = [ref_seq[i : i + 3] for i in range(0, length, 3)]
    alt_codons = [alt_seq[i : i + 3] for i in range(0, length, 3)]

    S, N = gene_sites(ref_seq[:length])
    Sd, Nd = 0.0, 0.0
    for c1, c2 in zip(ref_codons, alt_codons):
        if CODON_TABLE.get(c1, "*") == "*" or CODON_TABLE.get(c2, "*") == "*":
            continue
        sd, nd = codon_pathway_diff(c1, c2)
        Sd += sd
        Nd += nd

    ps = Sd / S if S > 0 else float("nan")
    pn = Nd / N if N > 0 else float("nan")
    return pn, ps, (pn / ps if ps > 0 else float("nan")), S, N, Sd, Nd


def trna_lookup():
    rows = list(csv.DictReader(open(RESULTS / "near_cognate_map.tsv"), delimiter="\t"))
    return {r["codon"]: int(r["final_tRNA_count_for_coAI"]) for r in rows}


def main():
    trna = trna_lookup()
    manifest = list(csv.DictReader(open(BASE / "manifest.tsv"), delimiter="\t"))
    patients = sorted(set(r["patient_num"] for r in manifest), key=int)
    successors = [p for p in patients if p != "1"]

    predecessor_seqs = {}
    for seg in ["L", "M", "S"]:
        path = BASE / "orfs" / "by_patient" / f"Patient_1_{seg}.fasta"
        if path.exists():
            _, seq = read_fasta(path)[0]
            predecessor_seqs[seg] = strip_stop(seq)

    pnps_rows = []
    mutation_rows = []

    for p in successors:
        for seg in ["L", "M", "S"]:
            path = BASE / "orfs" / "by_patient" / f"Patient_{p}_{seg}.fasta"
            if not path.exists() or seg not in predecessor_seqs:
                continue
            header, alt_seq = read_fasta(path)[0]
            alt_seq = strip_stop(alt_seq)
            ref_seq = predecessor_seqs[seg]
            length_mismatch = len(ref_seq) != len(alt_seq)

            pn, ps, ratio, S, N, Sd, Nd = pn_ps(ref_seq, alt_seq)
            pnps_rows.append(
                {
                    "patient_num": p,
                    "segment": seg,
                    "length_mismatch_note": "TRUNCATED_ALT_SEQ_COMPARED_ON_SHARED_LENGTH" if length_mismatch else "",
                    "S_sites": round(S, 3),
                    "N_sites": round(N, 3),
                    "Sd": round(Sd, 3),
                    "Nd": round(Nd, 3),
                    "ps": round(ps, 6) if ps == ps else "NA",
                    "pn": round(pn, 6) if pn == pn else "NA",
                    "pn_ps_ratio": round(ratio, 6) if ratio == ratio else "NA",
                }
            )

            # per-codon mutation detail
            length = min(len(ref_seq), len(alt_seq))
            length -= length % 3
            for i in range(0, length, 3):
                c1, c2 = ref_seq[i : i + 3], alt_seq[i : i + 3]
                if c1 == c2:
                    continue
                aa1, aa2 = CODON_TABLE.get(c1), CODON_TABLE.get(c2)
                if aa1 is None or aa2 is None:
                    continue
                syn = aa1 == aa2
                pol1 = POLARITY.get(aa1, "NA")
                pol2 = POLARITY.get(aa2, "NA")
                trna1 = trna.get(c1, "NA")
                trna2 = trna.get(c2, "NA")
                mutation_rows.append(
                    {
                        "patient_num": p,
                        "segment": seg,
                        "position_nt": i + 1,
                        "codon_ref": c1,
                        "codon_alt": c2,
                        "aa_ref": aa1,
                        "aa_alt": aa2,
                        "synonymous": syn,
                        "polarity_ref": pol1,
                        "polarity_alt": pol2,
                        "polarity_change": "" if syn else (f"{pol1}->{pol2}" if pol1 != pol2 else "none"),
                        "tRNA_count_ref_codon": trna1,
                        "tRNA_count_alt_codon": trna2,
                        "tRNA_count_change": (trna2 - trna1) if isinstance(trna1, int) and isinstance(trna2, int) else "NA",
                    }
                )

    with open(RESULTS / "pnps_per_patient_segment.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pnps_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(pnps_rows)
    print(f"Wrote {len(pnps_rows)} pn/ps rows (successors x segments)")

    if mutation_rows:
        with open(RESULTS / "mutations_detail.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(mutation_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(mutation_rows)
        print(f"Wrote {len(mutation_rows)} individual codon-change rows")
    else:
        print("No codon changes detected anywhere - writing empty placeholder")
        (RESULTS / "mutations_detail.tsv").write_text(
            "patient_num\tsegment\tposition_nt\tcodon_ref\tcodon_alt\taa_ref\taa_alt\tsynonymous\tpolarity_ref\tpolarity_alt\tpolarity_change\ttRNA_count_ref_codon\ttRNA_count_alt_codon\ttRNA_count_change\n"
        )


if __name__ == "__main__":
    main()
