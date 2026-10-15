#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""R3.1/R3.2 diagnostic — are the branch-model omega contrasts tracking
divergence depth rather than biology?

Motivation: the M0 F61 fit gives a tree length of 475 substitutions per codon
and one internal branch pinned at codeml's upper bound t = 50. At that depth dS
is saturated, so omega on a long branch is pushed toward zero for purely
numerical reasons, while a short, shallow branch retains a recoverable dS and
therefore a larger apparent omega. If that is happening, then EVERY branch
partition that happens to isolate shallow branches will look like it has
elevated omega -- which would explain, without any biology, (i) why the New
World HPS clade (the shallowest clade) is the best-fitting single partition,
(ii) why the pathogenic terminal branches look elevated under F61, and
(iii) why F3x4 and F61 disagree about the sign.

Test: per-species omega (from the 22 single-foreground codeml runs) against
that species' own terminal branch length, and against its total root-to-tip
path length. Spearman, because the relationship need not be linear.
"""
import json
import numpy as np
from scipy import stats
from Bio import Phylo

A = ("/PATH/TO/"
     "Project_Hantavirus/Reviews_ans/analysis/R3.1_deconfound_2026-10-12")

PATH = {"HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV",
        "CHOV", "LANV"}
RODENT = PATH | {"PHV", "TULV"}
NW = {"ANDV", "BAYV", "BCCV", "CHOV", "LANV", "SNV"}

d = json.load(open(f"{A}/results/04_fits.json"))
pertip = d["pertip"]

tree = Phylo.read(f"{A}/inputs/S22_ml_clean.nwk", "newick")
tiplen, roottip = {}, {}
for t in tree.get_terminals():
    tiplen[t.name] = t.branch_length
    roottip[t.name] = tree.distance(tree.root, t)

out = []
def say(s=""):
    out.append(s)
    print(s)

say("# Is the branch-level signal divergence depth rather than biology?")
say()
say("Per-species omega from the 22 single-foreground codeml runs, against that")
say("species' terminal branch length and root-to-tip distance on the ML tree.")
say()
say("| species | category | host | NW clade | terminal branch | root-to-tip | omega F3x4 | omega F61 |")
say("|---|---|---|---|---|---|---|---|")
rows = []
for sp in sorted(tiplen, key=lambda s: tiplen[s]):
    w3 = pertip["F3x4"].get(sp)
    w6 = pertip["F61"].get(sp)
    rows.append((sp, tiplen[sp], roottip[sp], w3, w6))
    say(f"| {sp} | {'pathogenic' if sp in PATH else 'environmental'} | "
        f"{'rodent' if sp in RODENT else 'non-rodent'} | "
        f"{'yes' if sp in NW else ''} | {tiplen[sp]:.4f} | {roottip[sp]:.4f} | "
        f"{'%.5f' % w3 if w3 is not None else 'NA'} | "
        f"{'%.5f' % w6 if w6 is not None else 'NA'} |")
say()

say("## Correlations")
say()
say("| omega from | vs | Spearman rho | p | n |")
say("|---|---|---|---|---|")
res = {}
for cf in ("F3x4", "F61"):
    for lab, dd in (("terminal branch length", tiplen),
                    ("root-to-tip distance", roottip)):
        pairs = [(dd[sp], pertip[cf][sp]) for sp in dd
                 if pertip[cf].get(sp) is not None]
        x = np.array([p[0] for p in pairs])
        y = np.array([p[1] for p in pairs])
        rho, p = stats.spearmanr(x, y)
        res[f"{cf}|{lab}"] = dict(rho=float(rho), p=float(p), n=len(x))
        say(f"| {cf} | {lab} | {rho:+.3f} | {p:.4g} | {len(x)} |")
say()

say("## Branch lengths by group (ML tree, substitutions per codon)")
say()
say("| group | n | mean terminal branch | median |")
say("|---|---|---|---|")
for nm, members in (("pathogenic", PATH),
                    ("environmental", set(tiplen) - PATH),
                    ("New World HPS clade", NW),
                    ("rodent-borne", RODENT),
                    ("non-rodent", set(tiplen) - RODENT)):
    v = np.array([tiplen[s] for s in members])
    say(f"| {nm} | {len(v)} | {v.mean():.4f} | {np.median(v):.4f} |")
say()
a = np.array([tiplen[s] for s in PATH])
b = np.array([tiplen[s] for s in set(tiplen) - PATH])
U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
say(f"Pathogenic vs environmental terminal branch length: Mann-Whitney "
    f"U = {U:.1f}, p = {p:.4g}.")
a = np.array([tiplen[s] for s in NW])
b = np.array([tiplen[s] for s in set(tiplen) - NW])
U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
say(f"New World HPS clade vs the rest: Mann-Whitney U = {U:.1f}, p = {p:.4g}.")
say()
with open(f"{A}/results/09_depth_artifact.json", "w") as fh:
    json.dump(res, fh, indent=1)
with open(f"{A}/results/09_depth_artifact.md", "w") as fh:
    fh.write("\n".join(out) + "\n")
print("\nwrote results/09_depth_artifact.{md,json}")
