#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""R3.1 — step 1.

(a) Independently verify the topology facts this analysis leans on:
      - is any single bipartition == the 10 pathogenic species?
      - is any single bipartition == the 12 rodent-borne species?
      - what is the smallest clade containing all 10 pathogenic species, and
        the largest clade that is purely pathogenic?
(b) Emit codeml branch-labelled trees for the head-to-head branch models.

All labellings mark TERMINAL branches only (ancestral host/disease states are
unknown), which is the same convention as the R2.2 two-ratio run, so the
models are directly comparable.
"""
import itertools
import sys
from io import StringIO

from Bio import Phylo

A = "/PATH/TO/Project_Hantavirus/Reviews_ans/analysis/R3.1_deconfound_2026-10-12"

PATHOGENIC = {"HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV",
              "CHOV", "LANV"}
RODENT = PATHOGENIC | {"PHV", "TULV"}
NW_HPS = {"ANDV", "BAYV", "BCCV", "CHOV", "LANV", "SNV"}

tree = Phylo.read(f"{A}/inputs/S22_ml_clean.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
assert len(tips) == 22, tips
assert PATHOGENIC <= tips and RODENT <= tips

out = []
def say(s=""):
    out.append(s)
    print(s)

say("# R3.1 step 1 — topology verification (independent re-check)")
say()
say(f"tips: {len(tips)}")
say(f"pathogenic (n={len(PATHOGENIC)}): {' '.join(sorted(PATHOGENIC))}")
say(f"rodent-borne (n={len(RODENT)}): {' '.join(sorted(RODENT))}")
say()

# --- bipartitions of the unrooted tree -------------------------------------
# For every internal branch, the set of tips on one side is one half of a
# bipartition. A set S is "separated by a single bipartition" iff S or its
# complement equals one of those halves.
clades = []
for cl in tree.get_nonterminals():
    s = frozenset(t.name for t in cl.get_terminals())
    clades.append(s)
bips = set()
for s in clades:
    bips.add(s)
    bips.add(frozenset(tips - s))

def separated(S):
    return frozenset(S) in bips

say("## Is the group separated by a single bipartition?")
say(f"- pathogenic monophyletic (10 spp):   {separated(PATHOGENIC)}")
say(f"- rodent-borne monophyletic (12 spp): {separated(RODENT)}")
say(f"- New World HPS clade (6 spp):        {separated(NW_HPS)}")
say()

# --- smallest clade containing all pathogenic ------------------------------
cands = [s for s in clades if PATHOGENIC <= s]
if cands:
    sm = min(cands, key=len)
    say(f"## Smallest clade containing ALL 10 pathogenic species: {len(sm)} tips")
    say(f"   {' '.join(sorted(sm))}")
    say(f"   extra (non-pathogenic) species dragged in: "
        f"{' '.join(sorted(sm - PATHOGENIC))}")
else:
    say("## No clade in the rooted tree contains all 10 pathogenic species")
say()

# largest purely-pathogenic clade
pure = [s for s in clades if s <= PATHOGENIC]
if pure:
    lg = max(pure, key=len)
    say(f"## Largest purely pathogenic clade: {len(lg)}/10 tips")
    say(f"   {' '.join(sorted(lg))}")
say()

# smallest clade containing all rodent-borne
cands = [s for s in clades if RODENT <= s]
if cands:
    sm = min(cands, key=len)
    say(f"## Smallest clade containing ALL 12 rodent-borne species: {len(sm)} tips")
    say(f"   {' '.join(sorted(sm))}")
    say(f"   non-rodent species dragged in: {' '.join(sorted(sm - RODENT))}")
say()

say("## All internal clades (rooted as given), size >= 2")
for s in sorted(clades, key=lambda x: (len(x), sorted(x))):
    tag = []
    if s <= PATHOGENIC:
        tag.append("all-pathogenic")
    if s <= RODENT:
        tag.append("all-rodent")
    if not (s & PATHOGENIC):
        tag.append("no-pathogenic")
    say(f"  n={len(s):2d}  {' '.join(sorted(s)):<60s} {','.join(tag)}")

with open(f"{A}/results/01_topology_check.md", "w") as fh:
    fh.write("\n".join(out) + "\n")

# --- labelled trees for codeml ---------------------------------------------
raw = open(f"{A}/inputs/S22_ml_clean.nwk").read().strip()

def label_tips(nwk, mapping):
    """Insert ' #k' after each named tip according to mapping {tip: k}."""
    s = nwk
    for tip, k in mapping.items():
        if k == 0:
            continue
        assert f"{tip}:" in s, tip
        s = s.replace(f"{tip}:", f"{tip} #{k}:", 1)
    return s

trees = {}
# H_disease : pathogenic tips = #1 (reproduces R2.2's S22_pathtips.nwk)
trees["disease2"] = label_tips(raw, {t: 1 for t in PATHOGENIC})
# H_host    : rodent-borne tips = #1  (10 pathogenic + PHV + TULV)
trees["host2"] = label_tips(raw, {t: 1 for t in RODENT})
# H_3class  : pathogenic tips #1, rodent-borne environmental tips (PHV,TULV) #2
m = {t: 1 for t in PATHOGENIC}
m.update({"PHV": 2, "TULV": 2})
trees["three"] = label_tips(raw, m)
# R3.3c sensitivity: TULV counted as pathogenic (11 foreground tips)
trees["disease2_tulvpath"] = label_tips(raw, {t: 1 for t in PATHOGENIC | {"TULV"}})
# R3.3c sensitivity: TULV dropped from the disease labelling entirely is the
# same as H_disease on a tree where TULV is its own class
m = {t: 1 for t in PATHOGENIC}
m["TULV"] = 2
trees["disease2_tulvown"] = label_tips(raw, m)
# clade alternative: the New World HPS clade incl. its internal branches.
# codeml propagates a label to a clade when placed after the closing paren.
nw_clade_sub = "(((SNV #1:0.2604487865,((ANDV #1:0.1848953969,LANV #1:0.1878804831) #1:0.0509631718,CHOV #1:0.2215584809) #1:0.0541645618) #1:0.0447237692,(BCCV #1:0.2144142898,BAYV #1:0.1267293429) #1:0.0933002608)"
orig_sub = "(((SNV:0.2604487865,((ANDV:0.1848953969,LANV:0.1878804831):0.0509631718,CHOV:0.2215584809):0.0541645618):0.0447237692,(BCCV:0.2144142898,BAYV:0.1267293429):0.0933002608)"
assert orig_sub in raw, "NW clade substring not found — tree changed"
trees["cladeNW"] = raw.replace(orig_sub, nw_clade_sub, 1)

for name, nwk in trees.items():
    with open(f"{A}/inputs/tree_{name}.nwk", "w") as fh:
        fh.write(nwk + "\n")
    print(f"wrote inputs/tree_{name}.nwk  (#1 count={nwk.count('#1')}, "
          f"#2 count={nwk.count('#2')})")

# sanity: disease2 must be identical in labelling to R2.2's S22_pathtips.nwk
ref = open("/PATH/TO/"
           "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12/"
           "inputs/S22_pathtips.nwk").read().strip()
print("disease2 identical to R2.2 S22_pathtips.nwk:", trees["disease2"] == ref)
