# `Revision_2026/` — follow-up analysis code

Everything here was written after the original analysis was complete. It is
kept separate from `Epuyen_outbreak/` and `Genus_Diversity/`, which are the
original pipeline, so that the two can be told apart: if you want to know what
the original pipeline did, read those; if you want to know how a reported number
was produced or checked, read this.

Three corrections were made *inside* the original pipeline rather than here,
because they change what that pipeline computes. They are listed in section 9 of
the top-level README.

---

## Read this before running anything

**Absolute paths were redacted.** These scripts were written to run in one
private working directory, with inputs laid out as they were on that machine.
Every absolute path has been rewritten to `/PATH/TO/`, and each affected file
says so at the top. The code is otherwise unchanged from the version that
produced the reported values — nothing was tidied, renamed or re-run for
publication.

**The inputs are not here.** This repository holds code only. The alignments,
trees, PAML outputs and intermediate tables these scripts read are in the data
deposit (see the DOI in the top-level README). Several scripts are also
read-only consumers of outputs from the original pipeline — they say so in
their own docstrings, and they do not write back into it.

**External tools.** The scripts in `R2.2_phylo_dnds_*`, `R3.1_deconfound_*` and
`R4.4_Msegment_*` drive IQ-TREE 3.1.1 and PAML 4.10.10 (`yn00` and `codeml`),
which must be on `PATH`. The `.ctl` files beside them are the control files as
used, and are part of the method record: they carry the settings the Methods
section states (`cleandata = 0`, `NSsites = 0`, `model = 2`, kappa estimated,
starting dN/dS 0.4).

---

## What each directory is

| directory | what it does |
|---|---|
| `CAI_metTrp_recompute_2026-10-14/` | The CAI correction. Recomputes every CAI value with the single-codon Met and Trp families excluded (Sharp & Li 1987), and refuses to emit anything until the earlier values are first reproduced exactly with those families included. Also harmonises the lung panel's B2M record with the kidney panel's, emits the corrected Table 7 and Supplementary Tables S1 and S7, and identifies — by reproduction, from the earlier p-values alone — what that version's "adjustment for GC3" actually was. |
| `R2.2_phylo_dnds_2026-10-12/` | The maximum-likelihood tree and the PAML work. Prepares the 22-species nucleoprotein codon alignment, infers the ML tree, then runs `yn00` (NG86 and YN00) and `codeml` branch models on it. Includes the codon-frequency sweep and the three-category branch model. This is where the corrected divergence estimators come from; the original pipeline's `02_pairwise_dnds.py` computes uncorrected pn/ps only. |
| `AUDIT_numbers_2026-10-12/` | An adversarial audit of the directory above. Every parser is written from scratch against the raw PAML output rather than reusing the original parser, so that a mistake in the original would surface as a disagreement instead of being reproduced. Covers the site sets, the saturation counts, branch support, topology, the free-ratio model and the sweep. |
| `R2.2_interpretation_2026-10-12/` | Whether the pathogenic-versus-environmental difference is explained by divergence rather than by selection. Read-only on the pairwise table. Also re-parses the free-ratio model with a boundary correction. |
| `R3.1_deconfound_2026-10-12/` | Rival explanations for the same contrast, on the same 231 pairs: relabelling by reservoir host class instead of by disease, by New World clade, by study effort, by record vintage and by sequencing depth, plus likelihood-ratio tests, a phylogenetic ANOVA and convergence checks. |
| `R3.1_origins_on_ML_tree_2026-10-12/` | Counts independent origins of human pathogenicity on the 22-species nucleoprotein tree. Written for a second paper that was abandoned; the result was folded into this one. |
| `R4.4_Msegment_dnds_2026-10-12/` | The M segment (glycoprotein precursor), added later. Mirrors the nucleoprotein analysis so the two are directly comparable, with the same alignment diagnostics, the same NG86-versus-YN00 comparison, the same saturation diagnostics and the same rival partitions. |
| `PINPS_within_species_2026-10-12/` | Within-species nucleoprotein pn/ps, as an independent check at a level where saturation is not a problem. Builds and purity-checks a per-species sequence inventory from GenBank, aligns, computes pn/ps, validates the NG86 implementation, and assesses power. The largest set of scripts here, because most of the work is inventory construction. |
| `R2.3a_reservoir_proxy_cai_2026-10-12/` | Viral CAI measured against a reservoir-host codon usage reference rather than a human one. The named reservoir species has no annotated genome, so *Peromyscus maniculatus* is used as an explicit proxy. |
| `R2.3b_isolate_host_cai_2026-10-12/` | Viral CAI of human-derived against animal-derived isolates, tested within each viral species, with host class harvested from the isolate records. |
| `R2.3c_tissue_cai_2026-10-12/` | CAI weight tables derived from the genes most expressed in a tissue, rather than from a generic human table — lung, kidney, and 49 further tissues as controls. This is the analysis that failed its own positive controls, which is why the paper withdraws its tissue-tropism inferences. |
| `R4.5_R4.8_outbreak_counts_coai_2026-10-12/` | Verification for the outbreak section: exact comparison counts, the tRNA co-adaptation index, a nucleotide-identity and ambiguity-handling fix, a tree-support audit and the codon-frequency sweep summary. |

## Negative results are kept

Four of these directories report that something could not be concluded:
`R2.3a` and `R2.3b` (host-reference CAI is dominated by the reference table's
GC3 content), `R2.3c` (no tissue resolution, by its own positive controls) and
parts of `R3.1_deconfound` (reservoir host class fits the data at least as well
as disease association). They are here for the same reason they are in the
paper: the attempt was made, and reporting it as uninformative is not the same
as not having made it.
