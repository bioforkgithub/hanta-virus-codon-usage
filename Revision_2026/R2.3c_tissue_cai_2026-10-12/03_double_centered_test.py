#!/usr/bin/env python3
"""R2.3c + R4.1 step 3 - the properly controlled form of the question.

Step 2 showed every query set (virus, housekeeping panels, tissue-ENRICHED
positive controls) produces the SAME tissue ranking: retina top, lung ~40th,
kidney ~18th of 51. That means raw per-tissue CAI is dominated by a property of
the REFERENCE (how concentrated that tissue's top-expressed genes are on
optimal codons), which moves every query up or down together, and carries no
query-specific information.

The signal R2.3c actually asks about is therefore the INTERACTION: after
removing the tissue main effect (some references give everyone high CAI) and
the sequence main effect (some sequences are optimal everywhere), is the
residual for lung/kidney positive for the virus?

  resid[i][t] = cai[i][t] - rowmean[i] - colmean[t] + grandmean

The positive control decides whether the residual has any resolution:
lung-ENRICHED human genes MUST peak at lung, kidney-ENRICHED at kidney. If
they do and the virus does not, the negative viral result is real. If they do
not, no CAI-based tissue-tropism statistic is recoverable at all and that is
the answer to both R2.3c and R4.1.
"""
import collections
import csv
import glob
import json
import math
import os
import random
import re
import sys

from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cai_lib  # noqa: E402

HPA, CDSFA, W1, RESULTS, LABELS = sys.argv[1:6]
TOPN = 500
random.seed(20261012)

gene_re = re.compile(r"gene:(ENSG\d+)")
bio_re = re.compile(r"transcript_biotype:(\S+)")
sym_re = re.compile(r"gene_symbol:(\S+)")
best = {}


def consider(hdr, seq):
    g, b = gene_re.search(hdr), bio_re.search(hdr)
    if not g or not b or b.group(1) != "protein_coding":
        return
    seq = seq.upper()
    if len(seq) < 300 or len(seq) % 3 or set(seq) - set("ACGT"):
        return
    if not seq.startswith("ATG") or seq[-3:] not in cai_lib.STOPS:
        return
    cur = best.get(g.group(1))
    if cur is None or len(seq) > len(cur[1]):
        s = sym_re.search(hdr)
        best[g.group(1)] = (s.group(1) if s else "", seq)


hdr, buf = None, []
with open(CDSFA) as fh:
    for line in fh:
        if line.startswith(">"):
            if hdr:
                consider(hdr, "".join(buf))
            hdr, buf = line, []
        else:
            buf.append(line.strip())
if hdr:
    consider(hdr, "".join(buf))
gcounts = {g: cai_lib.codon_counts(s) for g, (_s, s) in best.items()}

expr = collections.defaultdict(dict)
with open(HPA) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    i_g, i_t, i_v = head.index("Gene"), head.index("Tissue"), head.index("nTPM")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) > max(i_g, i_t, i_v):
            try:
                expr[p[i_t]][p[i_g]] = float(p[i_v])
            except ValueError:
                pass
tissues = sorted(expr)

W = {}
for t in tissues:
    r = sorted(((v, g) for g, v in expr[t].items() if g in gcounts), reverse=True)
    tot = collections.Counter()
    for _v, g in r[:TOPN]:
        tot.update(gcounts[g])
    W[t] = cai_lib.weights_proper(cai_lib.fractions_from_counts(tot))

# ------------------------------------------------------------ sequence sets
labels = {}
with open(LABELS) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        labels[row["accession"]] = row

N_PAT = re.compile(r"nucleocapsid|nucleoprotein", re.I)
GPC_PAT = re.compile(r"glycoprotein|precursor|\bGPC\b", re.I)
L_PAT = re.compile(r"polymerase|\bL protein\b", re.I)
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
                    sets["genus22_N"]["%s_%s" % (labels[rec.id]["abbrev"], rec.id)] = s

for p in sorted(glob.glob(os.path.join(W1, "Epuyen_outbreak", "gbk", "*"))):
    for rec in SeqIO.parse(p, "genbank"):
        for f in rec.features:
            prod = (f.qualifiers.get("product") or [""])[0] or f.type
            s = clean(f.extract(rec.seq))
            if not s:
                continue
            if N_PAT.search(prod):
                sets["epuyen_ANDV_N"][rec.id] = s
            elif L_PAT.search(prod):
                sets["epuyen_ANDV_L"][rec.id] = s
            elif GPC_PAT.search(prod):
                sets["epuyen_ANDV_GPC"][rec.id] = s

for nm, path in (("lung_HKG_panel", os.path.join(W1, "Lung_HKG", "lung_hkg_cds.fasta")),
                 ("kidney_HKG_panel", os.path.join(W1, "Kidney_HKG", "kidney_hkg_cds.fasta"))):
    if os.path.exists(path):
        for rec in SeqIO.parse(path, "fasta"):
            s = clean(rec.seq)
            if s:
                sets[nm][rec.id] = s

CTRL_TISSUES = ("lung", "kidney", "liver", "testis", "skeletal muscle", "pancreas",
                "cerebral cortex", "spleen")
for t in CTRL_TISSUES:
    out = []
    for g in gcounts:
        v = expr[t].get(g, 0.0)
        if v < 10:
            continue
        others = [expr[o].get(g, 0.0) for o in tissues if o != t]
        m = sum(others) / len(others)
        if v >= 5 * max(m, 1e-9) and v == max([v] + others):
            out.append((v, g))
    out.sort(reverse=True)
    sets["ENRICHED_" + re.sub(r"\W+", "_", t)] = {
        (best[g][0] or g): best[g][1] for g in [g for _v, g in out[:100]]}

# ------------------------------------------------------------------- matrix
rows, rowset = [], []
for name, seqs in sorted(sets.items()):
    for sid, s in sorted(seqs.items()):
        vals = []
        ok = True
        for t in tissues:
            c, _u, _k = cai_lib.cai_proper(s, W[t])
            if c is None:
                ok = False
                break
            vals.append(c)
        if ok:
            rows.append(vals)
            rowset.append((name, sid))

nr, nc = len(rows), len(tissues)
grand = sum(sum(r) for r in rows) / (nr * nc)
rowmean = [sum(r) / nc for r in rows]
colmean = [sum(rows[i][j] for i in range(nr)) / nr for j in range(nc)]
resid = [[rows[i][j] - rowmean[i] - colmean[j] + grand for j in range(nc)]
         for i in range(nr)]

with open(os.path.join(RESULTS, "double_centered_residuals.tsv"), "w") as fh:
    fh.write("set\tsequence\t" + "\t".join(tissues) + "\n")
    for i, (nm, sid) in enumerate(rowset):
        fh.write("%s\t%s\t" % (nm, sid) + "\t".join("%.6f" % v for v in resid[i]) + "\n")


def group_profile(name):
    idx = [i for i, (nm, _s) in enumerate(rowset) if nm == name]
    prof = [sum(resid[i][j] for i in idx) / len(idx) for j in range(nc)]
    order = sorted(range(nc), key=lambda j: -prof[j])
    rank = {tissues[order[k]]: k + 1 for k in range(nc)}
    # bootstrap over sequences in the group
    reps, lr, kr, home = 400, [], [], []
    hometis = None
    if name.startswith("ENRICHED_"):
        hometis = name[len("ENRICHED_"):].replace("_", " ")
        if hometis not in tissues:
            hometis = None
    for _ in range(reps):
        bs = [random.choice(idx) for _ in idx]
        p = [sum(resid[i][j] for i in bs) / len(bs) for j in range(nc)]
        o = sorted(range(nc), key=lambda j: -p[j])
        names = [tissues[k] for k in o]
        lr.append(names.index("lung") + 1)
        kr.append(names.index("kidney") + 1)
        if hometis:
            home.append(names.index(hometis) + 1)
    lr.sort()
    kr.sort()
    home.sort()
    res = {
        "n": len(idx),
        "lung_resid": prof[tissues.index("lung")],
        "kidney_resid": prof[tissues.index("kidney")],
        "lung_rank": rank["lung"], "kidney_rank": rank["kidney"],
        "lung_rank_ci95": [lr[int(.025 * reps)], lr[int(.975 * reps)]],
        "kidney_rank_ci95": [kr[int(.025 * reps)], kr[int(.975 * reps)]],
        "top3": [tissues[order[k]] for k in range(3)],
        "bottom3": [tissues[order[nc - 1 - k]] for k in range(3)],
        "max_abs_resid": max(abs(v) for v in prof),
    }
    if hometis:
        res["home_tissue"] = hometis
        res["home_rank"] = rank[hometis]
        res["home_rank_ci95"] = [home[int(.025 * reps)], home[int(.975 * reps)]]
        res["home_resid"] = prof[tissues.index(hometis)]
    return res


profiles = {nm: group_profile(nm) for nm in sorted(sets)}
out = {
    "n_rows": nr, "n_tissues": nc, "topN_reference_genes": TOPN,
    "grand_mean_cai": grand,
    "tissue_main_effect_range": max(colmean) - min(colmean),
    "sequence_main_effect_range": max(rowmean) - min(rowmean),
    "residual_sd": math.sqrt(sum(v * v for r in resid for v in r) / (nr * nc)),
    "profiles": profiles,
}
with open(os.path.join(RESULTS, "R2.3c_double_centered_summary.json"), "w") as fh:
    json.dump(out, fh, indent=2, sort_keys=True)

print("tissue main effect range (raw CAI): %.5f" % out["tissue_main_effect_range"])
print("sequence main effect range (raw CAI): %.5f" % out["sequence_main_effect_range"])
print("residual SD after double-centering: %.5f" % out["residual_sd"])
print()
print("%-28s %4s %9s %9s %10s %10s  %-22s %s" % (
    "group", "n", "lungRank", "kidRank", "lungResid", "kidResid", "home tissue (control)",
    "homeRank [95% CI]"))
for nm, p in profiles.items():
    home = ""
    if "home_rank" in p:
        home = "%-22s %d %s" % (p["home_tissue"], p["home_rank"], p["home_rank_ci95"])
    print("%-28s %4d %4d%-5s %4d%-5s %+10.5f %+10.5f  %s" % (
        nm, p["n"], p["lung_rank"], str(p["lung_rank_ci95"]),
        p["kidney_rank"], str(p["kidney_rank_ci95"]),
        p["lung_resid"], p["kidney_resid"], home))
print()
print("POSITIVE CONTROL verdict: does each ENRICHED_<tissue> set rank its own tissue first?")
for nm, p in profiles.items():
    if "home_rank" in p:
        print("  %-28s home=%-16s rank %d of %d  %s" % (
            nm, p["home_tissue"], p["home_rank"], nc,
            "PASS" if p["home_rank"] == 1 else "FAIL"))
