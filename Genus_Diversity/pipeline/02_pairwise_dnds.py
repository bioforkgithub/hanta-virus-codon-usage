"""
Pairwise pn/ps (Nei & Gojobori 1986 counting method) across all C(22,2)=231
species pairs of the across-species nucleoprotein codon alignment, plus per-
species CAI, ENC and GC3. Classifies each pair as pathogenic-pathogenic,
environmental-environmental, or cross (pathogenic-environmental) using the
manifest's category column.

LABELLING, corrected at revision. What this script computes is pn and ps as
uncorrected proportions and their ratio. No correction is applied for multiple
hits, for transition/transversion bias or for codon usage, so the quantity is
pn/ps and not dN/dS; the submitted version of the manuscript called it dN/dS,
which was inaccurate. The column and the output file are named accordingly
here (they were `dNdS` and `pairwise_dnds_matrix.tsv` in v1.1).

This matters in practice, not only in name: ps reaches or exceeds 0.75 - the
value at which the Jukes-Cantor correction becomes undefined - in 204 of the
231 pairs on the site sets this estimator uses. The corrected estimators
reported in the revised manuscript come from PAML and are in
`Revision_2026/phylogeny_and_paml/`.

Reuses the same synonymous/nonsynonymous site- and difference-counting
functions already validated on the Epuyen outbreak data
(Epuyen_outbreak/pipeline/06_pnps_mutations.py), applied here to codon
columns from the MAFFT protein-guided alignment with any gapped column
(in either sequence of the pair) excluded - the standard way to run NG86
on a pairwise basis from a multiple alignment.
"""
import csv
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "Epuyen_outbreak" / "pipeline"))
from shared_utils import read_fasta, load_cai_weights, compute_cai, compute_enc, strip_stop

EPUYEN_PIPELINE = Path(__file__).parent.parent.parent / "Epuyen_outbreak" / "pipeline"
sys.path.insert(0, str(EPUYEN_PIPELINE))
from importlib.util import spec_from_file_location, module_from_spec

spec = spec_from_file_location("pnps_mod", EPUYEN_PIPELINE / "06_pnps_mutations.py")
pnps_mod = module_from_spec(spec)
spec.loader.exec_module(pnps_mod)

BASE = Path(__file__).parent.parent
RESULTS = BASE / "results"
RESULTS.mkdir(exist_ok=True)


def gc3(seq):
    """Fraction of G or C at the third position, over sense codons only."""
    seq = strip_stop(seq)
    n = g = 0
    for i in range(0, len(seq) - 2, 3):
        c = seq[i : i + 3]
        aa = pnps_mod.CODON_TABLE.get(c)
        if aa is None or aa == "*":
            continue
        n += 1
        if c[2] in "GC":
            g += 1
    return g / n if n else float("nan")


def pairwise_ng86(codons_a, codons_b):
    """NG86 Sd/Nd/S/N over codon columns, skipping any column gapped in either sequence."""
    pairs = [(a, b) for a, b in zip(codons_a, codons_b) if "-" not in a and "-" not in b]
    S_total, N_total, Sd_total, Nd_total = 0.0, 0.0, 0.0, 0.0
    for c1, c2 in pairs:
        if pnps_mod.CODON_TABLE.get(c1, "*") == "*" or pnps_mod.CODON_TABLE.get(c2, "*") == "*":
            continue
        s1 = pnps_mod.synonymous_sites_per_codon(c1)
        s2 = pnps_mod.synonymous_sites_per_codon(c2)
        S_total += (s1 + s2) / 2
        N_total += 3 - (s1 + s2) / 2
        sd, nd = pnps_mod.codon_pathway_diff(c1, c2)
        Sd_total += sd
        Nd_total += nd
    ps = Sd_total / S_total if S_total > 0 else float("nan")
    pn = Nd_total / N_total if N_total > 0 else float("nan")
    ratio = pn / ps if ps > 0 else float("nan")
    return dict(n_codons_compared=len(pairs), S=S_total, N=N_total, Sd=Sd_total, Nd=Nd_total, ps=ps, pn=pn, pnps=ratio)


def main():
    manifest = list(csv.DictReader(open(BASE / "genus_diversity_manifest.tsv"), delimiter="\t"))
    category = {f"{r['abbrev']}|{r['species']}|{r['accession']}": r["category"] for r in manifest}
    abbrev_of = {f"{r['abbrev']}|{r['species']}|{r['accession']}": r["abbrev"] for r in manifest}

    aligned = read_fasta(BASE / "genus_diversity_codon_aligned.fasta")
    codon_seqs = {header: [seq[i : i + 3] for i in range(0, len(seq), 3)] for header, seq in aligned}
    headers = list(codon_seqs.keys())

    rows = []
    for h1, h2 in combinations(headers, 2):
        stats = pairwise_ng86(codon_seqs[h1], codon_seqs[h2])
        cat1, cat2 = category[h1], category[h2]
        comparison = "pathogenic-pathogenic" if cat1 == cat2 == "pathogenic" else (
            "environmental-environmental" if cat1 == cat2 == "environmental" else "cross"
        )
        rows.append({
            "species_A": abbrev_of[h1], "species_B": abbrev_of[h2],
            "category_A": cat1, "category_B": cat2, "comparison_type": comparison,
            "n_codons_compared": stats["n_codons_compared"],
            "S_sites": round(stats["S"], 3), "N_sites": round(stats["N"], 3),
            "Sd": round(stats["Sd"], 3), "Nd": round(stats["Nd"], 3),
            "ps": round(stats["ps"], 6) if stats["ps"] == stats["ps"] else "NA",
            "pn": round(stats["pn"], 6) if stats["pn"] == stats["pn"] else "NA",
            "pnps": round(stats["pnps"], 6) if stats["pnps"] == stats["pnps"] else "NA",
        })

    with open(RESULTS / "pairwise_pnps_matrix.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} pairwise pn/ps rows to results/pairwise_pnps_matrix.tsv")

    # --- per-species CAI/ENC/GC3 (own-species codon usage vs human high-expression reference) ---
    # CAI here excludes the single-codon Met and Trp families (Sharp & Li 1987);
    # see load_cai_weights. GC3 is written out because the revised manuscript
    # adjusts the CAI comparison between clinical categories for GC3 by
    # analysis of covariance - see 06_cai_gc3_ancova.R.
    weights = load_cai_weights()
    raw_recs = read_fasta(BASE / "genus_diversity_S_cds.fasta")
    cai_rows = []
    for header, seq in raw_recs:
        seq_clean = strip_stop(seq)
        cai = compute_cai(seq_clean, weights)
        enc = compute_enc(seq_clean)
        g3 = gc3(seq_clean)
        abbrev, species, acc = header.split("|")
        cai_rows.append({
            "abbrev": abbrev, "species": species, "accession": acc,
            "category": category[header], "length_nt": len(seq_clean),
            "CAI_vs_human_ref": round(cai, 4), "ENC": round(enc, 3),
            "GC3": round(g3, 6),
        })
    with open(RESULTS / "genus_diversity_cai_enc.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cai_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(cai_rows)
    print(f"Wrote {len(cai_rows)} per-species CAI/ENC/GC3 rows")


if __name__ == "__main__":
    main()
