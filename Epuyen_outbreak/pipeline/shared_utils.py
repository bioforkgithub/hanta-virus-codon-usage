"""
Shared utilities for the Epuyen Andes virus replication pipeline.

Implements CAI and ENC directly in Python rather than relying on this HPC's
EMBOSS install, because `cai` here is actually a symlink to codonW (a
different, incompatible program) and `chips` is a broken wrapper script
(verified 2026-07-30 - see WORK_PROGRESS.md). `cusp` was independently
confirmed genuine EMBOSS and is fine, but isn't needed once we compute
usage tables ourselves in Python.
"""
import csv
import random
from pathlib import Path
from math import log, exp

BASE = Path(__file__).parent.parent  # Epuyen_outbreak/

CODON_TABLE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "*", "TAG": "*",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "*", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}
SENSE_CODONS = [c for c in CODON_TABLE if CODON_TABLE[c] != "*"]
AA_TO_CODONS = {}
for c, aa in CODON_TABLE.items():
    if aa != "*":
        AA_TO_CODONS.setdefault(aa, []).append(c)

# The two single-codon amino-acid families, Met (ATG) and Trp (TGG). They have
# no synonyms, so they carry no information about codon choice; Sharp and Li
# (1987) exclude them from CAI. See `load_cai_weights`.
SINGLE_CODON_FAMILIES = sorted(aa for aa, cs in AA_TO_CODONS.items() if len(cs) == 1)
SINGLE_CODONS = {c for aa in SINGLE_CODON_FAMILIES for c in AA_TO_CODONS[aa]}

# amino-acid polarity classification (standard biochemistry)
POLARITY = {
    "Ala": "non-polar", "Val": "non-polar", "Leu": "non-polar", "Ile": "non-polar",
    "Pro": "non-polar", "Phe": "non-polar", "Trp": "non-polar", "Met": "non-polar",
    "Gly": "non-polar",
    "Ser": "polar", "Thr": "polar", "Cys": "polar", "Tyr": "polar", "Asn": "polar",
    "Gln": "polar",
    "Lys": "basic", "Arg": "basic", "His": "basic",
    "Asp": "acidic", "Glu": "acidic",
}


def read_fasta(path):
    """Return list of (header, sequence) tuples."""
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


def strip_stop(seq):
    seq = seq.upper().replace("U", "T")
    if len(seq) % 3 == 0 and seq[-3:] in ("TAA", "TAG", "TGA"):
        return seq[:-3]
    return seq


def load_cai_weights(cut_file=None, include_met_trp=False):
    """
    Load the Ehum.cut reference table (paper's cited Peden 1999 human
    high-expression codon usage set, bundled with EMBOSS) and compute
    proper CAI weights: w_ij = usage_ij / max(usage_i*) within each amino
    acid family (the normalization the original cai_fasta.py got wrong).

    `include_met_trp` selects whether the single-codon families Met and Trp
    enter the geometric mean:

      False (default) : they are omitted from the returned dictionary, so
                        `compute_cai` skips those codons entirely. This is the
                        standard definition of Sharp and Li (1987) and the
                        behaviour used for every CAI value in the revised
                        manuscript.
      True            : they are returned with w = 1, which is what the
                        submitted version of the manuscript did. A weight of
                        exactly one can only raise a geometric mean of terms
                        that are all <= 1, so those values were slightly
                        inflated. This mode is kept so the submitted numbers
                        can still be reproduced exactly from this repository,
                        which is how the correction was validated before any
                        value was changed.

    The correction is not a uniform offset: it depends on how much Met and Trp
    a sequence contains. The human reference panels fall by 0.0070-0.0081 CAI
    units, the Andes virus coding sequences by 0.0102-0.0113, and the 22
    nucleoprotein sequences by 0.0086-0.0175.
    """
    cut_file = cut_file or (BASE / "Ehum.cut")
    raw_freq = {}  # codon -> Frequency-per-1000 (column 3, "Number" is col 4 - we use Frequency to avoid rounding of Fraction)
    for line in open(cut_file):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        codon, aa, fraction, freq, number = parts[0], parts[1], parts[2], parts[3], parts[4]
        codon = codon.upper().replace("U", "T")
        if aa == "*":
            continue
        raw_freq[codon] = float(number)

    weights = {}
    for aa, codons in AA_TO_CODONS.items():
        if len(codons) == 1 and not include_met_trp:
            continue  # Met / Trp: no synonyms, excluded per Sharp and Li (1987)
        family_freqs = {c: raw_freq.get(c, 0.0) for c in codons}
        max_freq = max(family_freqs.values())
        for c, f in family_freqs.items():
            weights[c] = (f / max_freq) if max_freq > 0 else 1.0
    return weights


def compute_cai(seq, weights):
    seq = strip_stop(seq)
    if len(seq) % 3 != 0:
        raise ValueError("sequence length not a multiple of 3 after stop-strip")
    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    log_sum, n = 0.0, 0
    for c in codons:
        if c not in weights:
            # ambiguous / N-containing codons, and -- unless the weight table
            # was built with include_met_trp=True -- ATG and TGG
            continue
        w = weights[c]
        if w <= 0:
            w = 1e-6  # avoid log(0); shouldn't normally happen (own-codon w=1 always included)
        log_sum += log(w)
        n += 1
    if n == 0:
        raise ValueError("no valid codons found")
    return exp(log_sum / n)


def compute_enc(seq, undefined="as_published"):
    """
    Wright (1990) effective number of codons, "The 'effective number of
    codons' used in a gene", Gene 87:23-29 - same formula EMBOSS chips
    implements. Returns Nc in [20, 61].

    `undefined` selects what to do with a degeneracy class whose mean F cannot
    be estimated (no family in the class has n > 1) or estimates to exactly
    zero. F = 0 means a family whose codons were used perfectly evenly - with
    k codons used once each, sum(p^2) = 1/k = 1/n, so F = (n*sum(p^2) - 1)/(n - 1)
    is 0 - and 1/F is then a division by zero. It is a real, computed estimate,
    not a missing one, and it does arise here: the Ile family of the lung panel
    gene NDUFA3 has exactly one of each of ATT, ATC, ATA.

      'as_published' (default) : substitute F = 1 for such a class. This is what
            the submitted pipeline did, and it reproduces every ENC value
            printed in the manuscript, NDUFA3 (38.842) included. It is the
            default so that this repository reproduces the published tables.
            Be clear about what it means, though: F = 1 is MAXIMAL bias,
            imposed on a family whose codons are maximally even, so it pushes
            Nc down - the opposite direction from the truth. For NDUFA3 the
            affected class is the 3-fold one, contributing 1/F, so the error is
            bounded and small; it would be larger for the 2-fold (9/F) or
            4-fold (5/F) classes.
      'wright'                 : F3 = (F2 + F4)/2 when the 3-fold class is the
            unusable one - Wright's own suggestion for an absent or singleton
            Ile family. For NDUFA3 this gives 39.545 rather than 38.842.

    ENC is unaffected by the Met/Trp exclusion applied to CAI in the revision:
    single-codon families enter none of the degeneracy classes of this formula.
    """
    seq = strip_stop(seq)
    codons = [seq[i : i + 3] for i in range(0, len(seq), 3) if len(seq[i : i + 3]) == 3]
    counts = {}
    for c in codons:
        if c in CODON_TABLE and CODON_TABLE[c] != "*":
            counts[c] = counts.get(c, 0) + 1

    # group amino acids by degeneracy (fold) class; Ile is 3-fold, Met/Trp are 1-fold (excluded)
    fold_classes = {2: [], 3: [], 4: [], 6: []}
    for aa, aa_codons in AA_TO_CODONS.items():
        n = len(aa_codons)
        if n == 1:
            continue
        fold_classes.setdefault(n, []).append(aa_codons)

    def family_F(aa_codons):
        n_total = sum(counts.get(c, 0) for c in aa_codons)
        if n_total <= 1:
            return None  # insufficient data
        p = [counts.get(c, 0) / n_total for c in aa_codons]
        homo = sum(x * x for x in p)
        k = len(aa_codons)
        F = (n_total * homo - 1) / (n_total - 1)
        return F

    F_bar = {}
    for fold, families in fold_classes.items():
        vals = [family_F(fam) for fam in families]
        vals = [v for v in vals if v is not None]
        F_bar[fold] = sum(vals) / len(vals) if vals else None

    # Wright's formula: Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6.
    # Note the test below is `is None or == 0.0`, not a falsy test: a computed
    # F of exactly 0.0 is an estimate, and treating it as missing data was a
    # bug in the submitted version of this function. The two conventions are
    # documented in the docstring; 'as_published' keeps the published numbers.
    unusable = [k for k in (2, 3, 4, 6) if F_bar.get(k) is None or F_bar[k] == 0.0]
    if undefined == "as_published":
        for k in unusable:
            F_bar[k] = 1.0
    elif undefined == "wright":
        if 3 in unusable and F_bar.get(2) and F_bar.get(4):
            F_bar[3] = (F_bar[2] + F_bar[4]) / 2.0
            unusable.remove(3)
        for k in unusable:
            F_bar[k] = 1.0
    else:
        raise ValueError(f"unknown `undefined` convention: {undefined!r}")

    Nc = 2 + 9 / F_bar[2] + 1 / F_bar[3] + 5 / F_bar[4] + 3 / F_bar[6]
    return min(Nc, 61.0)


def windowed_mfe(seq, window=50, step=30, temp=37.0):
    """
    Paper's S1.2 method: mean minimum folding energy over sliding windows
    (50nt window, 30nt step, 37C), using seqfold (same tool the paper
    used - Timmons & Leshane 2019).
    """
    from seqfold import dg

    seq = seq.upper().replace("U", "T")
    energies = []
    i = 0
    n = len(seq)
    if n < window:
        energies.append(dg(seq, temp=temp))
    else:
        while i + window <= n:
            energies.append(dg(seq[i : i + window], temp=temp))
            i += step
        # trailing partial window (if any) to match paper's intent of
        # covering the whole CDS; paper doesn't specify a policy for the
        # remainder explicitly, so we include a final full-length window
        # anchored at the end to avoid dropping the tail entirely.
        if i < n and n - window >= 0 and (n - window) % step != 0:
            energies.append(dg(seq[n - window : n], temp=temp))
    return sum(energies) / len(energies)


def randomize_codons(seq, rng):
    """
    Paper's S1.1 method: synonymous codon randomization by sampling
    without replacement per amino acid within the gene. Preserves exact
    codon usage and protein sequence, only permutes which position gets
    which synonymous codon.
    """
    seq = seq.upper().replace("U", "T")
    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    aa_positions = {}
    for i, c in enumerate(codons):
        aa = CODON_TABLE.get(c, "?")
        aa_positions.setdefault(aa, []).append(i)

    new_codons = list(codons)
    for aa, positions in aa_positions.items():
        pool = [codons[i] for i in positions]
        rng.shuffle(pool)
        for i, c in zip(positions, pool):
            new_codons[i] = c
    return "".join(new_codons)
