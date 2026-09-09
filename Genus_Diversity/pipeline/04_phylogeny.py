"""
Neighbor-joining phylogeny of the 22 curated Orthohantavirus species from the
protein-guided alignment already built (01_translate_align.py), to visually
anchor the pathogenic-vs-environmental selective-constraint result (dN/dS,
02_pairwise_dnds.py) in an explicit evolutionary tree.

Uncorrected p-distance (proportion of differing aligned protein residues,
gap columns ignored per pair) is used to build the distance matrix; this is
a topology-illustration tree (to show the pathogenic/environmental split
qualitatively), not a publication-grade molecular-clock phylogeny.
"""
import csv
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).parent.parent
RESULTS = BASE / "results"
FIGURES = BASE / "figures"


def read_fasta(path):
    records = []
    header, seq = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header = line[1:]
            seq = []
        else:
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def p_distance(a, b):
    valid = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    if not valid:
        return 1.0
    diffs = sum(1 for x, y in valid if x != y)
    return diffs / len(valid)


def main():
    manifest = list(csv.DictReader(open(BASE / "genus_diversity_manifest.tsv"), delimiter="\t"))
    category = {r["abbrev"]: r["category"] for r in manifest}

    aligned = read_fasta(BASE / "genus_diversity_protein_aligned.fasta")
    abbrevs = [h.split("|")[0] for h, _ in aligned]
    seqs = [s for _, s in aligned]
    n = len(seqs)

    matrix = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            row.append(p_distance(seqs[i], seqs[j]) if i != j else 0.0)
        matrix.append(row)

    dm = DistanceMatrix(names=abbrevs, matrix=matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    tree.root_at_midpoint()

    Phylo.write(tree, str(RESULTS / "genus_diversity_nj_tree.nwk"), "newick")
    print(f"Wrote NJ tree ({n} taxa) to results/genus_diversity_nj_tree.nwk")

    # --- render, tip labels colored by clinical category ---
    fig, ax = plt.subplots(figsize=(8, 8))
    color_map = {"pathogenic": "#d64550", "environmental": "#2f7ed8"}

    def label_colors(name):
        if name is None:
            return "black"
        return color_map.get(category.get(str(name), ""), "black")

    Phylo.draw(
        tree, axes=ax, do_show=False,
        label_colors=label_colors,
        branch_labels=None,
        label_func=lambda c: c.name if c.is_terminal() else "",
    )
    ax.set_title("Neighbor-joining tree of 22 Orthohantavirus species (nucleoprotein)\n"
                  "red = pathogenic (HFRS/HPS-causing), blue = environmental (no known human disease)",
                  fontsize=11)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    plt.tight_layout()
    for ext, dpi in [("png", 150), ("pdf", 300), ("tiff", 300)]:
        fig.savefig(FIGURES / f"fig8_genus_phylogeny.{ext}", dpi=dpi)
    print("Wrote fig8_genus_phylogeny.{png,pdf,tiff}")


if __name__ == "__main__":
    main()
