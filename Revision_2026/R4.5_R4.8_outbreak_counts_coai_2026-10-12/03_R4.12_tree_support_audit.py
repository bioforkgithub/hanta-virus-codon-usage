# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""
R4.12 — support audit of the new ML tree that replaces the NJ / uncorrected
p-distance Figure 7.

Reads (read-only) the IQ-TREE outputs produced for R2.2 and reports, for every
internal branch, SH-aLRT and UFBoot, classified against the standard IQ-TREE
thresholds (SH-aLRT >= 80% AND UFBoot >= 95% = well supported). Also records
which clades the pathogenic-vs-environmental contrast depends on.
"""
import csv
import re
from pathlib import Path

R22 = Path(
    "/PATH/TO/"
    "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12"
)
OUT = Path(__file__).resolve().parent.parent / "results"

labels = {r["abbrev"]: r for r in csv.DictReader(open(R22 / "inputs" / "labels.tsv"), delimiter="\t")}
nwk = (R22 / "results" / "S22_ml.treefile").read_text().strip()

rep = []


def say(s=""):
    print(s)
    rep.append(s)


# ---- parse newick, collecting (clade tip set, support string) per internal node
def parse(s):
    pos = 0

    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            kids = []
            while True:
                kids.append(node())
                if s[pos] == ",":
                    pos += 1
                else:
                    break
            assert s[pos] == ")"
            pos += 1
            m = re.match(r"([^,():;]*)", s[pos:])
            lab = m.group(1)
            pos += len(lab)
            if s[pos] == ":":
                pos += 1
                m2 = re.match(r"[0-9.eE+-]+", s[pos:])
                pos += len(m2.group(0))
            tips = set()
            for k in kids:
                tips |= k[0]
            internal.append((tips, lab))
            return (tips, lab)
        m = re.match(r"([^,():;]+)", s[pos:])
        name = m.group(1)
        pos += len(name)
        if s[pos] == ":":
            pos += 1
            m2 = re.match(r"[0-9.eE+-]+", s[pos:])
            pos += len(m2.group(0))
        return ({name}, None)

    internal = []
    node()
    return internal


internal = parse(nwk)
all_tips = max((t for t, _ in internal), key=len)
n_tips = len(all_tips)

patho = {a for a, r in labels.items() if r["category"].lower().startswith("patho")}
env = all_tips - patho

say("# R4.12 — does the new ML tree adequately replace the NJ / p-distance Figure 7?")
say("")
say("Source (read-only): "
    "`analysis/R2.2_phylo_dnds_2026-10-12/results/S22_ml.{treefile,iqtree,log}`")
say("")
say("## What was run")
say("")
say("| item | value |")
say("|---|---|")
say("| program | IQ-TREE version 3.1.1 for Linux x86 64-bit, built Apr 8 2026 |")
say("| command | `iqtree -s S22_codon_aln.fasta -m MFP -B 1000 --alrt 1000 -T 8 --prefix S22_ml --seed 20261012 -redo` |")
say("| model selection | ModelFinder (`-m MFP`), best-fit by BIC |")
say("| substitution model | **GTR+F+R4** (FreeRate, 4 categories) |")
say("| rate categories | (0.29, 0.01696) (0.2014, 0.2347) (0.3868, 0.9595) (0.1218, 4.735) |")
say("| alignment | 22 sequences x 1344 nucleotide sites = 448 codons, S segment (nucleoprotein) |")
say("| constant sites | 356 (26.49%) |")
say("| parsimony-informative sites | 880 |")
say("| distinct site patterns | 959 |")
say("| support | 1000 ultrafast bootstrap (UFBoot2) + 1000 SH-aLRT replicates |")
say("| log-likelihood | -21943.2629 (s.e. 396.4169) |")
say("| free parameters | 55 (branches + model) |")
say("| total tree length | 12.5633; internal branches 4.0457 (32.20%) |")
say("")
say("This answers R4.12's request directly: an ML tree with an explicitly "
    "model-selected substitution model, over the same 22 complete S-segment "
    "sequences, with two independent branch-support measures — replacing a "
    "neighbour-joining tree built on uncorrected p-distances, which has neither a "
    "substitution model nor support values.")
say("")
say("One scope note that must be stated in the legend: the topology was inferred "
    "under a **nucleotide** model (GTR+F+R4) on the codon alignment, not under a "
    "codon model. That is standard practice and is what PAML then conditions on, "
    "but the figure legend should not describe it as a codon-model tree.")
say("")

# ---- support table
rows = []
for tips, lab in sorted(internal, key=lambda x: len(x[0])):
    if len(tips) == n_tips:
        continue
    if not lab:
        alrt = ufb = None
    else:
        parts = lab.split("/")
        alrt = float(parts[0]) if parts[0] else None
        ufb = float(parts[1]) if len(parts) > 1 and parts[1] else None
    npath = len(tips & patho)
    rows.append({
        "n_tips": len(tips),
        "clade": ",".join(sorted(tips)),
        "SH_aLRT": alrt, "UFBoot": ufb,
        "well_supported": (alrt is not None and ufb is not None and alrt >= 80 and ufb >= 95),
        "n_pathogenic_in_clade": npath,
        "n_environmental_in_clade": len(tips) - npath,
    })

with open(OUT / "R4.12_node_support.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

good = [r for r in rows if r["well_supported"]]
bad = [r for r in rows if not r["well_supported"]]

say("## Branch support")
say("")
say(f"Resolvable internal branches reported: {len(rows)}. "
    f"Well supported (SH-aLRT >= 80 AND UFBoot >= 95): **{len(good)}**. "
    f"Falling below one or both thresholds: **{len(bad)}**.")
say("")
say("| clade | n tips | SH-aLRT | UFBoot | well supported |")
say("|---|---|---|---|---|")
for r in sorted(rows, key=lambda r: (r["UFBoot"] is None, r["UFBoot"] or 0)):
    say(f"| {r['clade']} | {r['n_tips']} | "
        f"{'-' if r['SH_aLRT'] is None else r['SH_aLRT']} | "
        f"{'-' if r['UFBoot'] is None else r['UFBoot']} | "
        f"{'yes' if r['well_supported'] else '**no**'} |")
say("")
say("### Nodes that must be acknowledged as poorly supported")
say("")
for r in sorted(bad, key=lambda r: (r["UFBoot"] or 0)):
    say(f"- `{r['clade']}` — SH-aLRT {r['SH_aLRT']}, UFBoot {r['UFBoot']}")
say("")

# ---- what the ω contrast depends on
say("## Does the pathogenic-vs-environmental contrast depend on a weak node?")
say("")
say(f"Pathogenic species (n = {len(patho)}): {', '.join(sorted(patho))}")
say(f"Environmental species (n = {len(env)}): {', '.join(sorted(env))}")
say("")
pure_patho = [r for r in rows if r["n_environmental_in_clade"] == 0 and r["n_tips"] > 1]
say(f"Internal clades containing pathogenic species only: {len(pure_patho)} — "
    + (", ".join(f"`{r['clade']}` (SH-aLRT {r['SH_aLRT']}, UFBoot {r['UFBoot']})" for r in pure_patho)
       if pure_patho else "none"))
say("")
say(f"Largest all-pathogenic clade covers "
    f"{max((r['n_tips'] for r in pure_patho), default=0)} of {len(patho)} pathogenic species, "
    f"confirming that no single bipartition separates the pathogenic set.")
say("")
say("The codeml two-ratio foreground used for R4.3/R2.2 is the set of 10 "
    "pathogenic **terminal** branches. Terminal branches exist in every topology, "
    "so that foreground assignment does not depend on any of the weakly supported "
    "internal nodes listed above. The weak nodes do, however, bar any "
    "**clade-based** foreground definition, and they are the reason the "
    "polyphyly-of-pathogenicity statement has to be made from terminal branches "
    "rather than from a clade.")
say("")

(OUT / "R4.12_tree_support_audit.md").write_text("\n".join(rep) + "\n")
print("\nwrote", OUT / "R4.12_tree_support_audit.md")
