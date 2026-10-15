#!/usr/bin/env python3
"""Phase 4: how much of each species' 'polymorphism' is really between-lineage
divergence? Cluster at several identity thresholds and characterise the largest
cluster, which is the most defensible within-lineage sample for piN/piS.

Run: python3 -I 05_structure_sensitivity.py <inputs_dir> <results_dir>
Deterministic: no RNG; single-linkage, thresholds fixed below.
"""
import itertools
import json
import os
import sys
from collections import Counter

from Bio import SeqIO

ABBREVS = ["OXBV", "RKPV", "SWSV", "CBNV", "QZNV", "JJUV"]
THRESHOLDS = [0.90, 0.92, 0.95, 0.97]


def pdist(a, b):
    n = comp = 0
    for x, y in zip(a, b):
        if x in "ACGT" and y in "ACGT":
            comp += 1
            if x != y:
                n += 1
    return (n / comp) if comp else None


def clusters(accs, seqs, thresh):
    parent = {a: a for a in accs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in itertools.combinations(accs, 2):
        if len(seqs[a]) != len(seqs[b]):
            continue
        d = pdist(seqs[a], seqs[b])
        if d is not None and (1 - d) >= thresh:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    g = {}
    for a in accs:
        g.setdefault(find(a), []).append(a)
    return sorted((sorted(v) for v in g.values()), key=lambda c: (-len(c), c[0]))


def main():
    inputs_dir, results_dir = sys.argv[1], sys.argv[2]
    out = {}
    for ab in ABBREVS:
        fa = os.path.join(inputs_dir, ab, "N_cds_refined.fasta")
        seqs, info = {}, {}
        for r in SeqIO.parse(fa, "fasta"):
            acc = r.id.split("|")[0]
            seqs[acc] = str(r.seq).upper()
            # NB: use the FULL header (r.description), not r.id -- Biopython
            # truncates r.id at the first space, which silently cut host and
            # organism values such as 'Sorex araneus' down to 'Sorex'.
            parts = r.description.split("|")
            info[acc] = {"strain": parts[2] if len(parts) > 2 else "",
                         "host": next((p[5:] for p in parts if p.startswith("host=")), ""),
                         "country": next((p[8:] for p in parts if p.startswith("country=")), ""),
                         "org": next((p[4:] for p in parts if p.startswith("org=")), "")}
        accs = sorted(seqs)
        rec = {"n_independent": len(accs), "by_threshold": {}}

        for t in THRESHOLDS:
            cl = clusters(accs, seqs, t)
            big = cl[0] if cl else []
            dsb = []
            for a, b in itertools.combinations(big, 2):
                if len(seqs[a]) == len(seqs[b]):
                    d = pdist(seqs[a], seqs[b])
                    if d is not None:
                        dsb.append(d)
            dsb.sort()
            rec["by_threshold"][f"{t:.2f}"] = {
                "n_clusters": len(cl),
                "cluster_sizes": [len(c) for c in cl],
                "largest_n": len(big),
                "largest_n_haplotypes": len({seqs[a] for a in big}),
                "largest_members": big,
                "largest_organisms": dict(Counter(info[a]["org"] for a in big)),
                "largest_countries": dict(Counter(
                    info[a]["country"].split(":")[0].strip() for a in big)),
                "largest_hosts": dict(Counter(info[a]["host"] for a in big)),
                "largest_pdist_min": round(dsb[0], 5) if dsb else None,
                "largest_pdist_median": round(dsb[len(dsb) // 2], 5) if dsb else None,
                "largest_pdist_max": round(dsb[-1], 5) if dsb else None,
            }
        out[ab] = rec
        b95 = rec["by_threshold"]["0.95"]
        b90 = rec["by_threshold"]["0.90"]
        print(f"{ab}: n={len(accs)} | 0.90 -> {b90['n_clusters']} clusters, "
              f"largest {b90['largest_n']} (pdist max {b90['largest_pdist_max']}) | "
              f"0.95 -> {b95['n_clusters']} clusters, largest {b95['largest_n']} "
              f"(pdist max {b95['largest_pdist_max']})", flush=True)

    with open(os.path.join(results_dir, "structure_sensitivity_environmental-B.json"), "w") as fh:
        json.dump({"thresholds": THRESHOLDS, "species": out}, fh, indent=2, sort_keys=True)
    print("ok")


if __name__ == "__main__":
    main()
