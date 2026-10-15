"""Shared CAI machinery for R2.3b / R2.3c.

Two statistics are computed for every sequence, deliberately:

  cai_proper   Sharp & Li (1987) CAI. w(codon) = f(codon) / max f(synonym),
               so the optimal codon of every family has w = 1. Met, Trp and
               stop codons are EXCLUDED from the geometric mean (single-codon
               families carry no information).

  cai_ms       the statistic the original analysis actually computed
               (Work1/cai_fasta.py): the geometric mean of the raw
               per-amino-acid FRACTION, including Met and Trp at w = 1.
               Reproduced here only so the difference can be stated.

Nothing here reads the manuscript's results; weights are rebuilt from a codon
count or fraction table supplied by the caller.
"""
import math

from Bio.Data import CodonTable

TABLE = CodonTable.unambiguous_dna_by_id[1]
SYN = {}
for _c, _aa in TABLE.forward_table.items():
    SYN.setdefault(_aa, []).append(_c)
STOPS = set(TABLE.stop_codons)
SINGLETON_AA = {aa for aa, cs in SYN.items() if len(cs) == 1}  # M, W
CODON2AA = dict(TABLE.forward_table)


def read_fraction_table(path):
    """Read 'CODON fraction' lines (RNA or DNA alphabet) -> {DNA codon: float}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            out[parts[0].strip().upper().replace("U", "T")] = float(parts[1])
    return out


def fractions_from_counts(counts):
    """{codon: count} -> {codon: fraction within its amino-acid family}."""
    frac = {}
    for aa, cs in SYN.items():
        tot = sum(counts.get(c, 0) for c in cs)
        for c in cs:
            frac[c] = (counts.get(c, 0) / tot) if tot else 0.0
    return frac


def weights_proper(frac):
    """w = f / max f within family; informative families only."""
    w = {}
    for aa, cs in SYN.items():
        if aa in SINGLETON_AA:
            continue
        mx = max(frac.get(c, 0.0) for c in cs)
        if mx <= 0:
            continue
        for c in cs:
            w[c] = frac.get(c, 0.0) / mx
    return w


def codon_counts(seq):
    seq = seq.upper().replace("U", "T")
    counts = {}
    for i in range(0, len(seq) - len(seq) % 3, 3):
        c = seq[i:i + 3]
        if len(c) == 3 and not set(c) - set("ACGT"):
            counts[c] = counts.get(c, 0) + 1
    return counts


def _geomean(seq, w, skip_singletons, floor=0.01):
    logs, used, skipped = 0.0, 0, 0
    seq = seq.upper().replace("U", "T")
    for i in range(0, len(seq) - len(seq) % 3, 3):
        c = seq[i:i + 3]
        if len(c) != 3 or set(c) - set("ACGT"):
            skipped += 1
            continue
        if c in STOPS:
            skipped += 1
            continue
        aa = CODON2AA.get(c)
        if aa is None:
            skipped += 1
            continue
        if skip_singletons and aa in SINGLETON_AA:
            continue
        wc = w.get(c)
        if wc is None:
            skipped += 1
            continue
        logs += math.log(max(wc, floor))  # standard zero-weight floor
        used += 1
    if not used:
        return None, used, skipped
    return math.exp(logs / used), used, skipped


def cai_proper(seq, w):
    return _geomean(seq, w, skip_singletons=True)


def cai_manuscript_style(seq, frac):
    """Manuscript statistic: geometric mean of raw fractions, M/W included."""
    return _geomean(seq, frac, skip_singletons=False)
