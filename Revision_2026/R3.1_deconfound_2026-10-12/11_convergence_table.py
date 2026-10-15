#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""R3.1 — tabulate the convergence re-fits against the original runs.

A model whose lnL and omegas are reproduced from a different starting omega is
at a stable optimum and may be reported. A model that is not reproduced is at a
local optimum and its result must be withdrawn.
"""
import os
import re
import json

A = ("/PATH/TO/"
     "Project_Hantavirus/Reviews_ans/analysis/R3.1_deconfound_2026-10-12")
R22 = ("/PATH/TO/"
       "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12")


def parse(path):
    if not os.path.exists(path):
        return None
    txt = open(path).read()
    m = re.search(r"^lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*(-?[\d.]+)", txt, re.M)
    if not m:
        return None
    d = {"lnL": float(m.group(1))}
    w = re.search(r"^w \(dN/dS\) for branches:\s*(.+)$", txt, re.M)
    if w:
        d["w"] = [float(x) for x in w.group(1).split()]
    else:
        o = re.search(r"^omega \(dN/dS\) =\s*([\d.]+)", txt, re.M)
        d["w"] = [float(o.group(1))] if o else []
    return d


# model -> (original out.txt, [(start omega, conv out.txt), ...])
MODELS = {
    "M0 F61": (f"{R22}/results/codeml/m0_F61/out.txt",
               [("0.05", "m0_F61_w0.05"), ("1.0", "m0_F61_w1.0")]),
    "disease2 F61": (f"{R22}/results/codeml/two_ratio_F61/out.txt",
                     [("0.05", "disease2_F61_w0.05"),
                      ("1.0", "disease2_F61_w1.0")]),
    "host2 F61": (f"{A}/results/codeml/host2_F61/out.txt",
                  [("0.05", "host2_F61_w0.05"), ("1.0", "host2_F61_w1.0")]),
    "three F61": (f"{A}/results/codeml/three_F61/out.txt",
                  [("0.05", "three_F61_w0.05"), ("1.0", "three_F61_w1.0")]),
    "three_phvonly F61": (f"{A}/results/codeml/three_phvonly_F61/out.txt",
                          [("0.05", "three_phvonly_F61_w0.05")]),
    "cladeNW F61": (f"{A}/results/codeml/cladeNW_F61/out.txt",
                    [("0.05", "cladeNW_F61_w0.05")]),
    "disease2 F3x4": (f"{R22}/results/codeml/two_ratio_F3x4/out.txt",
                      [("1.0", "disease2_F3x4_w1.0")]),
    "three F3x4": (f"{A}/results/codeml/three_F3x4/out.txt",
                   [("1.0", "three_F3x4_w1.0")]),
}

out = []
def say(s=""):
    out.append(s)
    print(s)

say("# R3.1 — convergence check on the codeml fits")
say()
say("Every model below was originally fitted from starting omega 0.4. Each was")
say("re-fitted from the starting values shown. `dlnL` is the re-fit minus the")
say("original: a positive value means the original run was stuck at a WORSE")
say("optimum and the original result must not be reported.")
say()
say("| model | start omega | lnL | dlnL vs original | omegas | verdict |")
say("|---|---|---|---|---|---|")
verdict = {}
for name, (orig_path, runs) in MODELS.items():
    o = parse(orig_path)
    if o is None:
        say(f"| {name} | original | **MISSING** | | | |")
        continue
    say(f"| {name} | 0.4 (original) | {o['lnL']:.6f} | — | "
        f"{' '.join('%.5f' % x for x in o['w'])} | reference |")
    worst = 0.0
    anymissing = False
    for start, sub in runs:
        r = parse(f"{A}/results/codeml_conv/{sub}/out.txt")
        if r is None:
            say(f"| {name} | {start} | **PENDING** | | | |")
            anymissing = True
            continue
        d = r["lnL"] - o["lnL"]
        worst = max(worst, d)
        v = ("same optimum" if abs(d) < 0.01 else
             ("**BETTER optimum found — original was stuck**" if d > 0.01
              else "worse (original retained)"))
        say(f"| {name} | {start} | {r['lnL']:.6f} | {d:+.6f} | "
            f"{' '.join('%.5f' % x for x in r['w'])} | {v} |")
    verdict[name] = dict(max_gain=worst, pending=anymissing,
                         stable=bool(worst < 0.01 and not anymissing))
say()
say("## Verdict per model")
say()
say("| model | stable? | largest lnL gain from a different start |")
say("|---|---|---|")
for name, v in verdict.items():
    s = ("PENDING" if v["pending"] and v["max_gain"] < 0.01
         else ("yes" if v["stable"] else "**NO**"))
    say(f"| {name} | {s} | {v['max_gain']:+.6f} |")
say()
bad = [n for n, v in verdict.items() if not v["stable"] and not v["pending"]]
if bad:
    say(f"**Unstable fits: {', '.join(bad)}.** Any result that depends on these")
    say("must be withdrawn.")
else:
    say("No fit improved by more than 0.01 log-likelihood units from a")
    say("different starting omega, among the models checked so far.")
say()
with open(f"{A}/results/11_convergence.json", "w") as fh:
    json.dump(verdict, fh, indent=1)
with open(f"{A}/results/11_convergence.md", "w") as fh:
    fh.write("\n".join(out) + "\n")
print("\nwrote results/11_convergence.{md,json}")
