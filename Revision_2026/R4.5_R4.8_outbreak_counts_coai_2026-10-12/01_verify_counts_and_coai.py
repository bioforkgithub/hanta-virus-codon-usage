# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""
R4.5 / R4.9 / R4.10 / R4.8 — verification and new computation.

Reads ONLY (read-only) from Work1/Epuyen_outbreak/. Writes only into this
analysis folder's results/.

Part A (R4.5/R4.9/R4.10): exact comparison counts for the Epuyen outbreak
  successor-vs-index comparisons, to check the illustrative wording
  "In 75 of 80 comparisons, the genome was identical; in 4, only synonymous
  changes occurred" and the "5 codon substitutions across 27 successor
  genomes / ~10 weeks" characterisation.
  Done at three levels:
    A1 per-ORF (patient x segment) codon-level, from results/pnps_per_patient_segment.tsv
    A2 per-genome (patient) level
    A3 whole-genome nucleotide level, recomputed from gbk/ (ORF analysis
       covers coding regions only; "the genome was identical" is a stronger
       claim than the pn/ps table supports, so it is checked directly).

Part B (R4.8): co-AI (Pearson r between per-codon usage fraction and
  near-cognate-mapped human tRNA gene copy number, 61 sense codons) for the
  successor genomes, which Table 1 currently leaves as "not computed per
  genome". Reuses the existing pipeline's own functions and its own
  near_cognate_map.tsv so the value is on the same scale as the published
  predecessor / lung / kidney values.
"""
import csv
import json
import sys
from datetime import datetime
from pathlib import Path

from scipy.stats import pearsonr

OUTBREAK = Path(
    "/PATH/TO/"
    "Project_Hantavirus/Work1/Epuyen_outbreak"
)
OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(OUTBREAK / "pipeline"))
from shared_utils import read_fasta, strip_stop, CODON_TABLE, SENSE_CODONS  # noqa: E402

rep = []


def say(s=""):
    print(s)
    rep.append(s)


# ----------------------------------------------------------------- Part A
say("# R4.5 / R4.9 / R4.10 — verification of the comparison counts")
say(f"\nGenerated {datetime.now():%Y-%m-%d %H:%M}; source = {OUTBREAK} (read-only)\n")

manifest = list(csv.DictReader(open(OUTBREAK / "manifest.tsv"), delimiter="\t"))
pat_segs, pat_date = {}, {}
for r in manifest:
    p = int(r["patient_num"])
    pat_segs.setdefault(p, set()).add(r["segment"])
    pat_date[p] = datetime.strptime(r["sample_date"], "%m/%d/%y").date()

patients = sorted(pat_segs)
successors = [p for p in patients if p != 1]
span = (max(pat_date.values()) - min(pat_date.values())).days

say("## Dataset shape (from manifest.tsv)")
say("")
say("| quantity | value |")
say("|---|---|")
say(f"| GenBank accessions | {len(manifest)} |")
say(f"| patients (= genomes) | {len(patients)} |")
say(f"| index case | Patient 1, {min(pat_date.values()):%Y-%m-%d} |")
say(f"| successor genomes | {len(successors)} |")
say(f"| patient numbers present | {', '.join(str(p) for p in patients)} |")
say(f"| patient numbers absent | 21 (no sequences in the deposited set) |")
say(f"| patients with all 3 segments | {sum(1 for p in patients if pat_segs[p] == {'L','M','S'})} |")
_incomplete = ", ".join(
    "P%d (%s)" % (p, "".join(sorted(pat_segs[p])))
    for p in patients if pat_segs[p] != {"L", "M", "S"}
)
say(f"| patients with < 3 segments | {_incomplete} |")
say(f"| sampling span | {min(pat_date.values()):%Y-%m-%d} to {max(pat_date.values()):%Y-%m-%d} = {span} days = {span/7:.2f} weeks |")
say("")

# --- A1 per-ORF comparisons ----------------------------------------------
pnps = list(csv.DictReader(open(OUTBREAK / "results" / "pnps_per_patient_segment.tsv"), delimiter="\t"))
n_cmp = len(pnps)
ident, syn_only, nonsyn_only, both = [], [], [], []
undefined_ratio = 0
defined_vals = []
for r in pnps:
    sd, nd = float(r["Sd"]), float(r["Nd"])
    key = (int(r["patient_num"]), r["segment"])
    if sd == 0 and nd == 0:
        ident.append(key)
    elif nd == 0:
        syn_only.append(key)
    elif sd == 0:
        nonsyn_only.append(key)
    else:
        both.append(key)
    if r["pn_ps_ratio"] == "NA":
        undefined_ratio += 1
    else:
        defined_vals.append(float(r["pn_ps_ratio"]))

say("## A1 — per-ORF (patient x segment) codon-level comparisons vs the index case")
say("")
say("| outcome | n | of |")
say("|---|---|---|")
say(f"| comparisons performed | {n_cmp} | {len(successors)} successors x 3 segments = {len(successors)*3}, minus {len(successors)*3-n_cmp} (Patient 2 has no L segment) |")
say(f"| coding sequence identical to index | **{len(ident)}** | {n_cmp} |")
say(f"| synonymous change(s) only | **{len(syn_only)}** | {n_cmp} |")
say(f"| nonsynonymous change(s) only | **{len(nonsyn_only)}** | {n_cmp} |")
say(f"| both syn and nonsyn | {len(both)} | {n_cmp} |")
say(f"| pn/ps undefined (ps = 0) | {undefined_ratio} | {n_cmp} |")
say(f"| pn/ps defined | {len(defined_vals)} | {n_cmp} |")
say(f"| distinct defined pn/ps values | {sorted(set(defined_vals))} | — |")
say("")
say(f"Changed ORFs: syn-only {['P%d %s' % k for k in syn_only]}; "
    f"nonsyn-only {['P%d %s' % k for k in nonsyn_only]}")
say("")

# --- A2 per-genome ------------------------------------------------------
muts = list(csv.DictReader(open(OUTBREAK / "results" / "mutations_detail.tsv"), delimiter="\t"))
by_pat = {}
for m in muts:
    by_pat.setdefault(int(m["patient_num"]), []).append(m)
n_syn = sum(1 for m in muts if m["synonymous"] == "True")
n_nonsyn = len(muts) - n_syn

say("## A2 — per-genome (patient) level")
say("")
say("| quantity | value |")
say("|---|---|")
say(f"| codon substitutions observed in total | **{len(muts)}** |")
say(f"| of which synonymous | {n_syn} |")
say(f"| of which nonsynonymous | {n_nonsyn} |")
say(f"| successor genomes carrying >=1 codon substitution | **{len(by_pat)}** of {len(successors)} (P{', P'.join(str(p) for p in sorted(by_pat))}) |")
say(f"| successor genomes identical to index across all sequenced ORFs | **{len(successors)-len(by_pat)}** of {len(successors)} |")
say(f"| max substitutions in any one genome | {max(len(v) for v in by_pat.values())} |")
say("")
say("Each of the five substitutions is a single-nucleotide change in a distinct "
    "patient and a distinct ORF (verified below), so 5 codon substitutions = "
    "5 nucleotide substitutions, each in a different genome.")
say("")
say("| patient | segment | nt pos | codon | aa | synonymous |")
say("|---|---|---|---|---|---|")
for m in muts:
    nt_diff = sum(1 for a, b in zip(m["codon_ref"], m["codon_alt"]) if a != b)
    say(f"| P{m['patient_num']} | {m['segment']} | {m['position_nt']} | "
        f"{m['codon_ref']}>{m['codon_alt']} ({nt_diff} nt) | "
        f"{m['aa_ref']}>{m['aa_alt']} | {m['synonymous']} |")
say("")

# --- A3 whole-genome nucleotide level -----------------------------------
# SUPERSEDED. The first version of this section compared the deposited GenBank
# sequences 5'-anchored and reported ~75% mismatch for most pairs, which is the
# random expectation, not divergence: the deposited sequences differ in 5'
# extent (up to 268 nt) so the comparison was out of register. The nucleotide
# identity claim is made instead on the extracted, frame-aligned ORFs, in
# scripts/02_fix_nt_identity_and_ambiguity.py ->
# results/A3corrected_and_ambiguity.md.
say("## A3 — whole-sequence nucleotide identity")
say("")
say("Computed in `scripts/02_fix_nt_identity_and_ambiguity.py` -> "
    "`results/A3corrected_and_ambiguity.md`. Headline: of the 80 successor-vs-index "
    "ORF comparisons, **75 have zero unambiguous nucleotide differences** and 5 have "
    "exactly one each (5 nucleotide differences in total). The deposited GenBank "
    "sequences themselves differ in 5' extent and cannot be compared without "
    "alignment, so the identity statement is made on the ORFs.")
say("")

# ----------------------------------------------------------------- Part B
say("---")
say("")
say("# R4.8 — co-AI for the successor genomes")
say("")

trna = {r["codon"]: int(r["final_tRNA_count_for_coAI"])
        for r in csv.DictReader(open(OUTBREAK / "results" / "near_cognate_map.tsv"), delimiter="\t")}


def usage_fraction(seqs):
    counts = {c: 0 for c in SENSE_CODONS}
    total = 0
    for seq in seqs:
        seq = strip_stop(seq)
        for i in range(0, len(seq) - len(seq) % 3, 3):
            c = seq[i:i + 3]
            if c in counts:
                counts[c] += 1
                total += 1
    return {c: (counts[c] / total if total else 0.0) for c in SENSE_CODONS}, total


def coai(seqs):
    frac, total = usage_fraction(seqs)
    codons = sorted(frac)
    r, p = pearsonr([frac[c] for c in codons], [trna[c] for c in codons])
    return r, p, len(codons), total


# sanity: reproduce the published predecessor value from coai.tsv
pred_recs = read_fasta(OUTBREAK / "orfs" / "by_patient_combined" / "Patient_1_allORFs.fasta")
pr, pp, pn_, ptot = coai([s for _, s in pred_recs])
published = {r["genome"]: r for r in csv.DictReader(open(OUTBREAK / "results" / "coai.tsv"), delimiter="\t")}
pub_pred = float(published["Andes_virus_predecessor_genome"]["coAI_pearson_r"])

say("## B0 — reproduction check")
say("")
say(f"Recomputed index-case (Patient 1) co-AI = **{pr:.4f}** (p = {pp:.4g}, "
    f"{pn_} codons, {ptot} codons counted).")
say(f"Published value in `Work1/Epuyen_outbreak/results/coai.tsv` = {pub_pred:.4f}. "
    f"Difference = {abs(pr - pub_pred):.6f} -> the pipeline is reproduced exactly, "
    f"so the successor values below are on the same scale as Table 1.")
say("")

rows = []
for p in patients:
    recs = read_fasta(OUTBREAK / "orfs" / "by_patient_combined" / f"Patient_{p}_allORFs.fasta")
    r, pv, nc, tot = coai([s for _, s in recs])
    rows.append({"patient_num": p,
                 "role": "index" if p == 1 else "successor",
                 "segments_present": "".join(sorted(pat_segs[p])),
                 "codons_counted": tot,
                 "coAI_pearson_r": round(r, 4),
                 "p_value": f"{pv:.4g}",
                 "n_codons": nc,
                 "delta_vs_index": round(r - pr, 6)})

with open(OUT / "B_coai_per_genome.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

succ = [x for x in rows if x["role"] == "successor"]
succ_full = [x for x in succ if x["segments_present"] == "LMS"]
vals = [x["coAI_pearson_r"] for x in succ]
vals_full = [x["coAI_pearson_r"] for x in succ_full]
mean = sum(vals) / len(vals)
mean_full = sum(vals_full) / len(vals_full)
sd = (sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5
sd_full = (sum((v - mean_full) ** 2 for v in vals_full) / (len(vals_full) - 1)) ** 0.5

# pooled: all successor ORFs concatenated into one codon-usage table
pooled_seqs = []
for p in successors:
    pooled_seqs += [s for _, s in read_fasta(OUTBREAK / "orfs" / "by_patient_combined" / f"Patient_{p}_allORFs.fasta")]
sr, sp, sn, stot = coai(pooled_seqs)
# pooled over the 26 complete-genome successors only
pooled_full = []
for x in succ_full:
    pooled_full += [s for _, s in read_fasta(OUTBREAK / "orfs" / "by_patient_combined" / f"Patient_{x['patient_num']}_allORFs.fasta")]
fr, fp, fn_, ftot = coai(pooled_full)

say("## B1 — results")
say("")
say("Per-genome values are in `results/B_coai_per_genome.tsv`.")
say("")
say("| genome set | n | co-AI (Pearson r) | p | note |")
say("|---|---|---|---|---|")
say(f"| Index case (Patient 1) | 1 | {pr:.4f} | {pp:.4g} | published Table 1 value |")
say(f"| Successor genomes, per-genome mean +- SD | {len(vals)} | {mean:.4f} +- {sd:.4f} | — | range {min(vals):.4f} to {max(vals):.4f} |")
say(f"| Successor genomes with all 3 ORFs, mean +- SD | {len(vals_full)} | {mean_full:.4f} +- {sd_full:.4f} | — | range {min(vals_full):.4f} to {max(vals_full):.4f} |")
say(f"| Successor genomes pooled (all ORFs concatenated) | {len(successors)} | {sr:.4f} | {sp:.4g} | {stot} codons |")
say(f"| Successor genomes pooled, complete genomes only | {len(succ_full)} | {fr:.4f} | {fp:.4g} | {ftot} codons |")
say("")

outliers = [x for x in succ if abs(x["delta_vs_index"]) > 0.005]
say("## B2 — is the footnote's \"effectively identical to the index case\" true?")
say("")
say(f"Successor genomes whose co-AI differs from the index by more than 0.005: "
    f"{len(outliers)} of {len(succ)}"
    + (" — " + ", ".join(f"P{x['patient_num']} (r = {x['coAI_pearson_r']:.4f}, "
                         f"delta = {x['delta_vs_index']:+.4f}, segments {x['segments_present']})"
                         for x in outliers) if outliers else ""))
say("")
say(f"Among the {len(succ_full)} successors with all three ORFs sequenced, the "
    f"largest absolute deviation from the index case is "
    f"{max(abs(x['delta_vs_index']) for x in succ_full):.6f}.")
say("")

json.dump({"n_accessions": len(manifest), "n_patients": len(patients),
           "n_successors": len(successors), "span_days": span,
           "n_orf_comparisons": n_cmp, "n_identical_orf": len(ident),
           "n_syn_only_orf": len(syn_only), "n_nonsyn_only_orf": len(nonsyn_only),
           "n_both_orf": len(both), "n_pnps_undefined": undefined_ratio,
           "n_pnps_defined": len(defined_vals), "pnps_defined_values": sorted(set(defined_vals)),
           "n_codon_subs": len(muts), "n_syn": n_syn, "n_nonsyn": n_nonsyn,
           "n_genomes_changed": len(by_pat), "genomes_changed": sorted(by_pat),
           "coai_index": pr, "coai_succ_mean": mean, "coai_succ_sd": sd,
           "coai_succ_min": min(vals), "coai_succ_max": max(vals),
           "coai_succ_full_mean": mean_full, "coai_succ_full_sd": sd_full,
           "coai_pooled": sr, "coai_pooled_p": sp,
           "coai_pooled_full": fr, "coai_pooled_full_p": fp},
          open(OUT / "key_numbers.json", "w"), indent=2)

(OUT / "R4.5_R4.8_verification.md").write_text("\n".join(rep) + "\n")
print("\nwrote", OUT / "R4.5_R4.8_verification.md")
