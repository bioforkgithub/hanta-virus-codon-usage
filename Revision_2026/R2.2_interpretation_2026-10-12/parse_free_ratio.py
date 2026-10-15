#!/usr/bin/env python3
"""Read-only parse of the finished free-ratio codeml run; nothing is re-run."""
import json, re, statistics, sys
OUT, PHY, LAB, DEST = sys.argv[1:5]
# taxon order in the phylip file == codeml node numbering 1..n
order = []
with open(PHY) as fh:
    fh.readline()
    for line in fh:
        line = line.strip()
        if line and not set(line.upper()) <= set("ACGTN-? "):
            order.append(line.split()[0])
cat = {}
with open(LAB) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        cat[p[0]] = p[hdr.index("category")]
rows = []
started = False
for line in open(OUT):
    if re.match(r"\s*branch\s+t\s+N\s+S\s+dN/dS", line):
        started = True; continue
    if started:
        m = re.match(r"\s*(\d+)\.\.(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
        if not m:
            if rows: break
            continue
        a, b = int(m.group(1)), int(m.group(2))
        rows.append({"parent": a, "child": b, "t": float(m.group(3)),
                     "omega": float(m.group(6)), "dN": float(m.group(7)),
                     "dS": float(m.group(8)),
                     "taxon": order[b-1] if b <= len(order) else None})
for r in rows:
    r["category"] = cat.get(r["taxon"]) if r["taxon"] else "internal"
term = [r for r in rows if r["taxon"]]
path = [r for r in term if r["category"] == "pathogenic"]
envr = [r for r in term if r["category"] == "environmental"]
intl = [r for r in rows if not r["taxon"]]
def st(v, k="omega"):
    x = [r[k] for r in v]
    return {"n": len(x), "mean": statistics.mean(x), "median": statistics.median(x),
            "min": min(x), "max": max(x)}
res = {
  "n_branches": len(rows), "n_terminal": len(term), "n_internal": len(intl),
  "pathogenic_terminals": st(path), "environmental_terminals": st(envr),
  "internal_branches": st(intl),
  "pathogenic_terminals_dS": st(path, "dS"), "environmental_terminals_dS": st(envr, "dS"),
  "branches_omega_gt_1": sum(1 for r in rows if r["omega"] > 1),
  "branches_omega_gt_10": sum(1 for r in rows if r["omega"] > 10),
  "branches_dS_lt_0.01": sum(1 for r in rows if r["dS"] < 0.01),
  "branches_t_at_50_boundary": sum(1 for r in rows if r["t"] >= 49.999),
  "taxa_at_t_boundary": sorted(r["taxon"] for r in rows if r["t"] >= 49.999 and r["taxon"]),
  "branches": rows,
}
json.dump(res, open(DEST + "/free_ratio_branch_summary.json", "w"), indent=2, sort_keys=True)
for k, v in res.items():
    if k != "branches":
        print(k, "=", json.dumps(v))
