#!/usr/bin/env python3
"""R2.3b step 3 - viral CAI of human-derived vs animal-derived isolates.

Design: the contrast is tested WITHIN each viral species, because viral species
identity dominates codon usage and a pooled human-vs-animal test across species
would measure species composition, not host.

Reference codon usage: the same human table the manuscript used
(Work1/human_codon_usage_table.txt), so that the only thing that changes
relative to the earlier analysis is the grouping and the weight definition.
"""
import collections
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cai_lib  # noqa: E402

HUMAN_TABLE = sys.argv[1]
FASTA = sys.argv[2]
TSV = sys.argv[3]
OUTDIR = sys.argv[4]
MIN_N = 8  # per group, per species, to be reported as a test rather than a count

frac = cai_lib.read_fraction_table(HUMAN_TABLE)
w = cai_lib.weights_proper(frac)

meta = {}
with open(TSV) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        meta[row["accession"]] = row

seqs = {}
name = None
buf = []
with open(FASTA) as fh:
    for line in fh:
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
if name:
    seqs[name] = "".join(buf)

rows = []
for hdr, seq in seqs.items():
    acc = hdr.split("|")[0]
    m = meta.get(acc)
    if not m:
        continue
    cp, used_p, _ = cai_lib.cai_proper(seq, w)
    cm, used_m, _ = cai_lib.cai_manuscript_style(seq, frac)
    if cp is None:
        continue
    rows.append({
        "accession": acc, "organism": m["organism"], "group": m["group"],
        "host_evidence": m["host_evidence"], "country": m["country"],
        "coll_date": m["coll_date"], "strain": m["strain"],
        "cds_len": int(m["cds_len"]), "codons_used": used_p,
        "cai_proper": cp, "cai_manuscript_style": cm,
    })

with open(os.path.join(OUTDIR, "R2.3b_cai_per_sequence.tsv"), "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    wr.writeheader()
    for r in rows:
        wr.writerow(r)


# ---- statistics (exact Mann-Whitney U, normal approximation with tie
# correction; no scipy in this env, so the implementation is explicit) -------
def mannwhitney(a, b):
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return None
    allv = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks = [0.0] * (n1 + n2)
    i = 0
    ties = []
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = r
        ties.append(j - i + 1)
        i = j + 1
    r1 = sum(ranks[k] for k in range(len(allv)) if allv[k][1] == 0)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u2 = n1 * n2 - u1
    u = min(u1, u2)
    mu = n1 * n2 / 2.0
    n = n1 + n2
    tiesum = sum(t ** 3 - t for t in ties)
    sd = math.sqrt((n1 * n2 / 12.0) * ((n + 1) - tiesum / (n * (n - 1.0))))
    if sd == 0:
        return {"U": u, "p": 1.0, "z": 0.0, "A": u1 / (n1 * n2)}
    z = (u - mu) / sd
    p = math.erfc(abs(z) / math.sqrt(2))
    return {"U": u, "p": p, "z": z, "A": u1 / (n1 * n2)}  # A = P(human > animal)


def desc(v):
    v = sorted(v)
    n = len(v)
    mean = sum(v) / n
    med = v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2
    sd = math.sqrt(sum((x - mean) ** 2 for x in v) / (n - 1)) if n > 1 else 0.0
    return {"n": n, "mean": mean, "median": med, "sd": sd, "min": v[0], "max": v[-1]}


by_sp = collections.defaultdict(lambda: collections.defaultdict(list))
for r in rows:
    by_sp[r["organism"]][r["group"]].append(r["cai_proper"])

per_species = []
for sp, g in sorted(by_sp.items()):
    h, a = g.get("human", []), g.get("animal", [])
    rec = {"organism": sp, "n_human": len(h), "n_animal": len(a)}
    if len(h) >= MIN_N and len(a) >= MIN_N:
        rec["human"] = desc(h)
        rec["animal"] = desc(a)
        rec["test"] = mannwhitney(h, a)
        rec["delta_mean"] = rec["human"]["mean"] - rec["animal"]["mean"]
        rec["testable"] = True
    else:
        rec["testable"] = False
    per_species.append(rec)

testable = [r for r in per_species if r["testable"]]
# Stouffer combination over species with signed z (human - animal direction)
if testable:
    zs = []
    for r in testable:
        z = abs(r["test"]["z"]) * (1 if r["delta_mean"] > 0 else -1)
        zs.append(z)
    z_comb = sum(zs) / math.sqrt(len(zs))
    p_comb = math.erfc(abs(z_comb) / math.sqrt(2))
else:
    z_comb = p_comb = None

out = {
    "reference_table": HUMAN_TABLE,
    "n_sequences_scored": len(rows),
    "min_n_per_group_for_test": MIN_N,
    "per_species": per_species,
    "stouffer_across_testable_species": {
        "k_species": len(testable), "z": z_comb, "p": p_comb,
        "species": [r["organism"] for r in testable],
    },
    "weight_check": {
        "n_informative_codons": len(w),
        "max_weight_per_family_is_one": all(
            abs(max(w[c] for c in cs if c in w) - 1.0) < 1e-9
            for aa, cs in cai_lib.SYN.items() if aa not in cai_lib.SINGLETON_AA
        ),
    },
}
with open(os.path.join(OUTDIR, "R2.3b_stats.json"), "w") as fh:
    json.dump(out, fh, indent=2, sort_keys=True)

print(json.dumps(out["weight_check"], indent=2))
print("sequences scored:", len(rows))
print()
hdr = "%-34s %5s %5s | %-9s %-9s %-9s | %-9s %s"
print(hdr % ("organism", "nHum", "nAni", "meanHum", "meanAni", "delta", "p(MWU)", "A=P(H>A)"))
for r in per_species:
    if r["testable"]:
        print("%-34s %5d %5d | %9.5f %9.5f %+9.5f | %9.3g %.3f" % (
            r["organism"][:34], r["n_human"], r["n_animal"], r["human"]["mean"],
            r["animal"]["mean"], r["delta_mean"], r["test"]["p"], r["test"]["A"]))
    elif r["n_human"] + r["n_animal"] >= 5:
        print("%-34s %5d %5d | underpowered (needs >=%d per group)" % (
            r["organism"][:34], r["n_human"], r["n_animal"], MIN_N))
print()
print("Stouffer across %s testable species: z=%s p=%s" % (
    len(testable), z_comb, p_comb))
