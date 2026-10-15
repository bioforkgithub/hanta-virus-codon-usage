#!/usr/bin/env python3
"""
Break down WHY records were excluded, so "partials excluded" is a measured
number rather than a residual bucket.

Categories (mutually exclusive, in order):
  M_or_L_segment     : DEFINITION/features indicate M or L segment -> not our gene
  S_no_cds_feature   : S-segment (or N-gene) record with NO CDS annotation at all
                       -> recoverable by ORF-finding; counted separately
  N_cds_partial      : an N-like CDS exists but is too short / fuzzy / partial
  other_cds_only     : only non-N CDS (e.g. NSs alone)
  unclassifiable     : DEFINITION gives no segment

Usage: python3 -I diagnose_excluded.py <SPECIES_DIR>
"""
import sys, os, re, json
from Bio import SeqIO
from Bio.Seq import Seq

M_RE = re.compile(r"\bsegment M\b|\bM segment\b|glycoprotein|\bGPC\b|\bM gene\b", re.I)
L_RE = re.compile(r"\bsegment L\b|\bL segment\b|polymerase|\bRdRp\b|\bL gene\b", re.I)
S_RE = re.compile(r"\bsegment S\b|\bS segment\b|nucleocapsid|nucleoprotein|"
                  r"\bN gene\b|\bNP gene\b|\bS gene\b", re.I)
NSS_RE = re.compile(r"nonstructural|non-structural|\bNSs\b", re.I)

N_AA_MIN, N_AA_MAX = 380, 460


def longest_orf(seq):
    """Longest ATG..stop ORF on the given strand, returned as nt length."""
    s = str(seq).upper()
    best = 0
    for frame in range(3):
        i = frame
        while i < len(s) - 2:
            if s[i:i + 3] == "ATG":
                j = i
                while j < len(s) - 2:
                    c = s[j:j + 3]
                    if c in ("TAA", "TAG", "TGA"):
                        best = max(best, j + 3 - i)
                        break
                    j += 3
                else:
                    pass
            i += 3
    return best


def main():
    d = sys.argv[1]
    gb = os.path.join(d, "records.gb")
    tab = os.path.join(d, "records_table.tsv")
    # accessions already classified as having an N CDS candidate
    have_cds = set()
    for line in open(tab):
        p = line.rstrip("\n").split("\t")
        if p[0] == "accession":
            continue
        if p[3] not in ("0", ""):
            have_cds.add(p[0])

    cat = {}
    recoverable = []
    for rec in SeqIO.parse(gb, "genbank"):
        if rec.id in have_cds:
            continue
        defn = rec.description
        cdss = [f for f in rec.features if f.type == "CDS"]
        n_like = []
        for f in cdss:
            prod = " ".join(f.qualifiers.get("product", []))
            if NSS_RE.search(prod):
                continue
            try:
                L = len(f.extract(rec.seq))
            except Exception:
                continue
            n_like.append(L // 3)
        if M_RE.search(defn) and not S_RE.search(defn):
            c = "M_or_L_segment"
        elif L_RE.search(defn) and not S_RE.search(defn):
            c = "M_or_L_segment"
        elif S_RE.search(defn) or not defn.strip():
            if not cdss:
                orf = longest_orf(rec.seq)
                if N_AA_MIN * 3 <= orf <= N_AA_MAX * 3:
                    c = "S_no_cds_feature_FULL_ORF"
                    recoverable.append((rec.id, len(rec.seq), orf, defn))
                else:
                    c = "S_no_cds_feature_partial"
            elif n_like:
                c = "N_cds_partial"
            else:
                c = "other_cds_only"
        else:
            c = "unclassifiable"
        cat[c] = cat.get(c, 0) + 1

    out = dict(species_dir=os.path.basename(d.rstrip("/")),
               excluded_breakdown=cat,
               n_recoverable_full_orf=len(recoverable),
               recoverable_examples=sorted(recoverable)[:10])
    with open(os.path.join(d, "exclusion_breakdown.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
