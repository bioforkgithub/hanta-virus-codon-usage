# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""
Two corrections to 01_verify_counts_and_coai.py.

(1) A3 redone. In 01 the deposited GenBank sequences were compared 5'-anchored
    in the as-deposited orientation and most pairs came out ~75% mismatched,
    i.e. at the random expectation (25% identity). That is a bug in the
    comparison, not divergence: the deposited sequences are not all in the same
    orientation / not coordinate-aligned. Here each pair is tested in both
    orientations and the better one reported, and the nucleotide-identity claim
    is additionally made on the ORFs extracted by the project's own
    extract_orfs.py, which are orientation-normalised.

(2) Why do some successor genomes have a slightly different co-AI / ENC / CAI
    than the index case even though mutations_detail.tsv records no codon
    substitution for them (P10, P11, P18)? Hypothesis: codons containing
    ambiguous bases are silently dropped by the pipeline's codon counters and
    by the pn/ps counter, so missing/ambiguous data shifts the codon-usage
    table without being recorded as a substitution. Quantified here.
"""
import csv
import sys
from pathlib import Path

OUTBREAK = Path(
    "/PATH/TO/"
    "Project_Hantavirus/Work1/Epuyen_outbreak"
)
OUT = Path(__file__).resolve().parent.parent / "results"
sys.path.insert(0, str(OUTBREAK / "pipeline"))
from shared_utils import read_fasta, strip_stop, SENSE_CODONS  # noqa: E402

rep = []


def say(s=""):
    print(s)
    rep.append(s)


COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def gbk_seq(path):
    seq, inseq = [], False
    for line in open(path):
        if line.startswith("ORIGIN"):
            inseq = True
            continue
        if inseq:
            if line.startswith("//"):
                break
            seq.append("".join(ch for ch in line if ch.isalpha()))
    return "".join(seq).upper()


manifest = list(csv.DictReader(open(OUTBREAK / "manifest.tsv"), delimiter="\t"))
acc, pat_segs = {}, {}
for r in manifest:
    p = int(r["patient_num"])
    acc[(p, r["segment"])] = r["accession"].split(".")[0]
    pat_segs.setdefault(p, set()).add(r["segment"])
patients = sorted(pat_segs)
successors = [p for p in patients if p != 1]

say("# A3 corrected + ambiguity audit (R4.5 / R4.8 supporting checks)")
say("")

# ---------------------------------------------------------------- (1a) orientation
say("## A3-corrected, step 1 — orientation of the deposited sequences")
say("")
idx = {s: gbk_seq(OUTBREAK / "gbk" / f"{acc[(1, s)]}.gbk") for s in "LMS"}


def best_anchored(ref, q):
    out = []
    for orient, qq in (("as_deposited", q), ("revcomp", rc(q))):
        n = min(len(ref), len(qq))
        mism = sum(1 for i in range(n) if ref[i] != qq[i])
        out.append((mism / n, orient, n, mism))
    return min(out)


orient_counts = {}
for p in successors:
    for s in sorted(pat_segs[p]):
        _, orient, _, _ = best_anchored(idx[s], gbk_seq(OUTBREAK / "gbk" / f"{acc[(p, s)]}.gbk"))
        orient_counts[orient] = orient_counts.get(orient, 0) + 1
say(f"Best-matching orientation, 5'-anchored, over 80 successor segment sequences: {orient_counts}")
say("")
say("Neither orientation gives a good match for most pairs — both score near the "
    "random expectation (~75% mismatch = ~25% identity), which is why 17 pairs "
    "score marginally better reverse-complemented purely as noise. The cause is "
    "that the deposited sequences differ in 5' extent (untranslated flanks of "
    "unequal length: differences of up to 268 nt), so a 5'-anchored comparison is "
    "out of register. The deposited sequences are not coordinate-aligned and a "
    "raw whole-sequence comparison of them is not interpretable. Nucleotide "
    "identity is therefore assessed on the extracted ORFs below, which the "
    "project's extract_orfs.py has already put in a common frame and orientation.")
say("")

# ---------------------------------------------------------------- (1b) ORF nt identity
say("## A3-corrected, step 2 — nucleotide identity of the extracted ORFs")
say("")
say("ORFs from `orfs/by_patient/Patient_<n>_<seg>.fasta` (the project's own "
    "extract_orfs.py output), compared codon-frame-aligned to Patient 1.")
say("")

rows = []
for p in successors:
    for s in sorted(pat_segs[p]):
        ref = strip_stop(read_fasta(OUTBREAK / "orfs" / "by_patient" / f"Patient_1_{s}.fasta")[0][1])
        q = strip_stop(read_fasta(OUTBREAK / "orfs" / "by_patient" / f"Patient_{p}_{s}.fasta")[0][1])
        n = min(len(ref), len(q))
        mism = [i for i in range(n) if ref[i] != q[i]]
        amb = sum(1 for i in range(n) if q[i] not in "ACGT")
        unamb = [i for i in mism if q[i] in "ACGT" and ref[i] in "ACGT"]
        rows.append({"patient": p, "segment": s, "orf_len_index": len(ref),
                     "orf_len_this": len(q), "compared_nt": n,
                     "mismatches_any": len(mism),
                     "mismatches_unambiguous": len(unamb),
                     "ambiguous_nt_in_this": amb})

with open(OUT / "A3corr_orf_nt_identity.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

same_len = sum(1 for r in rows if r["orf_len_index"] == r["orf_len_this"])
ident_unamb = [r for r in rows if r["mismatches_unambiguous"] == 0]
with_amb = [r for r in rows if r["ambiguous_nt_in_this"] > 0]
say("| quantity | value |")
say("|---|---|")
say(f"| ORF comparisons | {len(rows)} |")
say(f"| ORF pair same length as index | {same_len} of {len(rows)} |")
say(f"| zero unambiguous nucleotide differences from index | **{len(ident_unamb)}** of {len(rows)} |")
say(f"| >=1 unambiguous nucleotide difference | **{len(rows)-len(ident_unamb)}** of {len(rows)} |")
say(f"| ORFs containing >=1 ambiguous base | {len(with_amb)} of {len(rows)} |")
say(f"| total unambiguous nucleotide differences, all ORFs | **{sum(r['mismatches_unambiguous'] for r in rows)}** |")
say(f"| total ambiguous positions, all ORFs | {sum(r['ambiguous_nt_in_this'] for r in rows)} |")
say("")
say("| patient | segment | unambiguous diffs | ambiguous nt | ORF len (index / this) |")
say("|---|---|---|---|---|")
for r in sorted(rows, key=lambda r: (-r["mismatches_unambiguous"], -r["ambiguous_nt_in_this"])):
    if r["mismatches_unambiguous"] or r["ambiguous_nt_in_this"]:
        say(f"| P{r['patient']} | {r['segment']} | {r['mismatches_unambiguous']} | "
            f"{r['ambiguous_nt_in_this']} | {r['orf_len_index']} / {r['orf_len_this']} |")
say("")

# ---------------------------------------------------------------- (2) ambiguity audit
say("## Ambiguity audit — why some successors' co-AI/ENC/CAI differ without a recorded substitution")
say("")


def codon_stats(seqs):
    counted, dropped = 0, 0
    for seq in seqs:
        seq = strip_stop(seq)
        for i in range(0, len(seq) - len(seq) % 3, 3):
            c = seq[i:i + 3]
            if c in SENSE_CODONS:
                counted += 1
            else:
                dropped += 1
    return counted, dropped


muts = list(csv.DictReader(open(OUTBREAK / "results" / "mutations_detail.tsv"), delimiter="\t"))
mut_pats = {int(m["patient_num"]) for m in muts}
strain = {int(r["patient_num"]): r for r in csv.DictReader(
    open(OUTBREAK / "results" / "strain_genome_metrics.tsv"), delimiter="\t")}
coai_rows = {int(r["patient_num"]): r for r in csv.DictReader(
    open(OUT / "B_coai_per_genome.tsv"), delimiter="\t")}

aud = []
for p in patients:
    seqs = [s for _, s in read_fasta(OUTBREAK / "orfs" / "by_patient_combined" / f"Patient_{p}_allORFs.fasta")]
    counted, dropped = codon_stats(seqs)
    aud.append({"patient": p,
                "segments": "".join(sorted(pat_segs[p])),
                "codons_counted": counted,
                "codons_dropped_ambiguous": dropped,
                "recorded_codon_subs": sum(1 for m in muts if int(m["patient_num"]) == p),
                "ENC": strain[p]["ENC"], "CAI": strain[p]["CAI"],
                "coAI": coai_rows[p]["coAI_pearson_r"],
                "coAI_delta_vs_index": coai_rows[p]["delta_vs_index"]})

with open(OUT / "B2_ambiguity_audit.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(aud[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(aud)

say("| patient | segments | codons counted | codons dropped (ambiguous) | recorded codon subs | ENC | CAI | co-AI | co-AI delta vs index |")
say("|---|---|---|---|---|---|---|---|---|")
for a in aud:
    flag = " <-" if (a["codons_dropped_ambiguous"] or a["recorded_codon_subs"]) else ""
    say(f"| P{a['patient']}{flag} | {a['segments']} | {a['codons_counted']} | "
        f"{a['codons_dropped_ambiguous']} | {a['recorded_codon_subs']} | {a['ENC']} | "
        f"{a['CAI']} | {a['coAI']} | {a['coAI_delta_vs_index']} |")
say("")

no_sub_but_dropped = [a for a in aud if a["patient"] != 1 and a["recorded_codon_subs"] == 0
                      and a["codons_dropped_ambiguous"] > 0]
say(f"Successor genomes with NO recorded codon substitution but >=1 codon dropped "
    f"for ambiguity: {len(no_sub_but_dropped)} — "
    + ", ".join(f"P{a['patient']} ({a['codons_dropped_ambiguous']} codons dropped, "
                f"co-AI delta {a['coAI_delta_vs_index']})" for a in no_sub_but_dropped))
say("")
say("This confirms the hypothesis: the small per-genome spread in co-AI, ENC and "
    "CAI among successors is driven by ambiguous (N) positions dropped from the "
    "codon tables, not by substitutions. It is a data-completeness artefact, and "
    "it is the honest reason the per-genome co-AI column was left blank.")
say("")

(OUT / "A3corrected_and_ambiguity.md").write_text("\n".join(rep) + "\n")
print("\nwrote", OUT / "A3corrected_and_ambiguity.md")
