#!/usr/bin/env python3
"""R2.3a - viral CAI under a human vs a reservoir-proxy codon usage reference.

The reservoir Oligoryzomys longicaudatus has 4 mRNA records
and no annotated genome in the genus, so it cannot supply a reference.
This script uses
*Peromyscus maniculatus* (GCF_003704035.1, the Sin Nombre virus reservoir and
the nearest cricetid with an annotated CDS set) as an EXPLICIT proxy.

Two matched reference pairs, so human and reservoir are built the same way:
  all_cds   codon usage over every clean CDS in the annotation
  ribosomal RPL*/RPS* genes, the standard stand-in for "highly expressed"
    when no expression atlas exists for the species

Reporting rule: a difference between two mammalian references is only
interpretable against how far apart the two references are in the first place,
so the reference-to-reference distance is reported alongside.
"""
import collections
import glob
import json
import math
import os
import re
import sys

from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cai_lib  # noqa: E402

PMAN_CDS, HUMAN_CDS, W1, LABELS, DEST = sys.argv[1:6]

RIBO = re.compile(r"\b(RPL|RPS)\d", re.I)


def load_counts(path, kind):
    """Return (all_cds_counts, ribosomal_counts, n_all, n_ribo)."""
    allc, ribc, na, nr = collections.Counter(), collections.Counter(), 0, 0
    for rec in SeqIO.parse(path, "fasta"):
        s = str(rec.seq).upper()
        if len(s) < 300 or len(s) % 3 or set(s) - set("ACGT"):
            continue
        if not s.startswith("ATG") or s[-3:] not in cai_lib.STOPS:
            continue
        c = cai_lib.codon_counts(s)
        allc.update(c)
        na += 1
        desc = rec.description
        if kind == "pman":
            m = re.search(r"\[gene=([^\]]+)\]", desc)
            gene = m.group(1) if m else ""
        else:
            m = re.search(r"gene_symbol:(\S+)", desc)
            gene = m.group(1) if m else ""
        if gene and RIBO.match(gene):
            ribc.update(c)
            nr += 1
    return allc, ribc, na, nr


pm_all, pm_rib, pm_na, pm_nr = load_counts(PMAN_CDS, "pman")
hs_all, hs_rib, hs_na, hs_nr = load_counts(HUMAN_CDS, "human")
sys.stderr.write("P.maniculatus CDS used: %d (ribosomal %d)\n" % (pm_na, pm_nr))
sys.stderr.write("H.sapiens CDS used: %d (ribosomal %d)\n" % (hs_na, hs_nr))

REFS = {
    "human_all_cds": cai_lib.weights_proper(cai_lib.fractions_from_counts(hs_all)),
    "human_ribosomal": cai_lib.weights_proper(cai_lib.fractions_from_counts(hs_rib)),
    "reservoirproxy_Pmaniculatus_all_cds":
        cai_lib.weights_proper(cai_lib.fractions_from_counts(pm_all)),
    "reservoirproxy_Pmaniculatus_ribosomal":
        cai_lib.weights_proper(cai_lib.fractions_from_counts(pm_rib)),
}

# ------------------------------------------------------------- viral sequences
labels = {}
with open(LABELS) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        labels[p[hdr.index("accession")]] = dict(zip(hdr, p))

N_PAT = re.compile(r"nucleocapsid|nucleoprotein", re.I)
sets = collections.defaultdict(dict)


def clean(s):
    s = str(s).upper()
    return None if (len(s) % 3 or set(s) - set("ACGT") or len(s) < 300) else s


for p in sorted(glob.glob(os.path.join(W1, "Genus_Diversity", "gbk", "*.gbk"))):
    for rec in SeqIO.parse(p, "genbank"):
        if rec.id not in labels:
            continue
        for f in rec.features:
            if f.type == "CDS" and N_PAT.search((f.qualifiers.get("product") or [""])[0]):
                s = clean(f.extract(rec.seq))
                if s:
                    key = "%s (%s)" % (labels[rec.id]["abbrev"], labels[rec.id]["category"])
                    sets["genus22_N"][key] = s
for p in sorted(glob.glob(os.path.join(W1, "Epuyen_outbreak", "gbk", "*"))):
    for rec in SeqIO.parse(p, "genbank"):
        for f in rec.features:
            prod = (f.qualifiers.get("product") or [""])[0] or f.type
            s = clean(f.extract(rec.seq))
            if s and N_PAT.search(prod):
                sets["epuyen_ANDV_N"][rec.id] = s

# ---------------------------------------------------------------------- score
def mean(v):
    return sum(v) / len(v)


rowsout = []
summary = {}
for setname, seqs in sorted(sets.items()):
    per_ref = {}
    for ref, w in REFS.items():
        vals = []
        for sid, s in seqs.items():
            c, _u, _k = cai_lib.cai_proper(s, w)
            if c is not None:
                vals.append((sid, c))
        per_ref[ref] = dict(vals)
    ids = sorted(set.intersection(*[set(d) for d in per_ref.values()]))
    for sid in ids:
        rowsout.append([setname, sid] + ["%.6f" % per_ref[r][sid] for r in sorted(REFS)])
    summary[setname] = {
        "n": len(ids),
        "mean_cai": {r: mean([per_ref[r][i] for i in ids]) for r in sorted(REFS)},
        "paired_delta_all_cds_human_minus_reservoir":
            mean([per_ref["human_all_cds"][i]
                  - per_ref["reservoirproxy_Pmaniculatus_all_cds"][i] for i in ids]),
        "paired_delta_ribosomal_human_minus_reservoir":
            mean([per_ref["human_ribosomal"][i]
                  - per_ref["reservoirproxy_Pmaniculatus_ribosomal"][i] for i in ids]),
        "n_sequences_higher_under_human_all_cds": sum(
            1 for i in ids if per_ref["human_all_cds"][i]
            > per_ref["reservoirproxy_Pmaniculatus_all_cds"][i]),
        "n_sequences_higher_under_human_ribosomal": sum(
            1 for i in ids if per_ref["human_ribosomal"][i]
            > per_ref["reservoirproxy_Pmaniculatus_ribosomal"][i]),
    }

with open(os.path.join(DEST, "R2.3a_cai_per_sequence.tsv"), "w") as fh:
    fh.write("set\tsequence\t" + "\t".join(sorted(REFS)) + "\n")
    for r in rowsout:
        fh.write("\t".join(r) + "\n")

# ------------------------------------------- how far apart are the references?
def rms(a, b):
    sh = [c for c in a if c in b]
    return math.sqrt(sum((a[c] - b[c]) ** 2 for c in sh) / len(sh))


refdist = {}
names = sorted(REFS)
for i, a in enumerate(names):
    for b in names[i + 1:]:
        refdist["%s|%s" % (a, b)] = rms(REFS[a], REFS[b])

out = {
    "reservoir_proxy": {
        "species": "Peromyscus maniculatus",
        "assembly": "GCF_003704035.1",
        "reason": "O. longicaudatus has 4 mRNA records and no annotated genome in genus",
        "n_cds_used": pm_na, "n_ribosomal_cds_used": pm_nr,
    },
    "human_reference": {"source": HUMAN_CDS, "n_cds_used": hs_na,
                        "n_ribosomal_cds_used": hs_nr},
    "sets": summary,
    "reference_to_reference_rms_weight_distance": refdist,
}
with open(os.path.join(DEST, "R2.3a_summary.json"), "w") as fh:
    json.dump(out, fh, indent=2, sort_keys=True)

print(json.dumps(out, indent=2, sort_keys=True))
