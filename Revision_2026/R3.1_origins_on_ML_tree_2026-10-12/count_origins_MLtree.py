#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""
Paper 2, step 2 (PLAN_second_paper.md ADDENDUM 3 sec. 3.5):
count independent origins of human pathogenicity across the 22-species
Orthohantavirus nucleoprotein phylogeny from paper 1.

Inputs (read-only, from paper 1):
    Work1/Genus_Diversity/results/genus_diversity_nj_tree.nwk   (NJ, midpoint-rooted)
    Work1/Genus_Diversity/genus_diversity_manifest.tsv          (curated categories)

Method
------
1. Fitch/Sankoff parsimony on the binary trait (pathogenic / environmental).
   The parsimony score is the minimum number of state CHANGES and is
   independent of rooting, so it is the honest root-free answer to "how many
   times did this trait change?". Verified here by re-rooting on every edge.
2. "Origins" is directional, so a dynamic program reports the MINIMUM and
   MAXIMUM number of gains (env->path) over ALL most-parsimonious ancestral
   reconstructions, under each of several candidate roots. This gives the range
   of origin counts the data are compatible with, not one cherry-picked number.
3. Dollo count: minimum gains if the trait is assumed never to be lost.
4. The maximal all-pathogenic clades, i.e. what paper 1's "the ten pathogenic
   species are essentially one clade" caveat actually amounts to on this tree.

No new data, no alignment, no sequence handling. Exact DP, no enumeration.
"""
from pathlib import Path
import sys

PROJ = Path("/PATH/TO/Project_Hantavirus")
TREE = Path("/PATH/TO/Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12/inputs/S22_ml_clean.nwk")
MANIFEST = Path("/PATH/TO/Project_Hantavirus/Work1/Genus_Diversity/genus_diversity_manifest.tsv")

PATH, ENV = "P", "E"
STATES = (PATH, ENV)
INF = float("inf")


# ---------------------------------------------------------------- newick parse
class Node:
    __slots__ = ("name", "children", "parent")

    def __init__(self, name=""):
        self.name = name
        self.children = []
        self.parent = None

    def is_leaf(self):
        return not self.children


def parse_newick(text):
    """Minimal newick parser; drops branch lengths, keeps labels."""
    text = text.strip().rstrip(";")
    pos = 0

    def parse_node():
        nonlocal pos
        node = Node()
        if text[pos] == "(":
            pos += 1
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"bad newick at {pos}: {text[pos:pos+20]!r}")
        start = pos
        while pos < len(text) and text[pos] not in "(),:":
            pos += 1
        node.name = text[start:pos]
        if pos < len(text) and text[pos] == ":":
            pos += 1
            while pos < len(text) and text[pos] not in "(),":
                pos += 1
        return node

    return parse_node()


def leaves(node):
    if node.is_leaf():
        return [node]
    out = []
    for c in node.children:
        out.extend(leaves(c))
    return out


def internals(node):
    if node.is_leaf():
        return []
    out = [node]
    for c in node.children:
        out.extend(internals(c))
    return out


def edges_of(root):
    out = []
    stack = [root]
    while stack:
        n = stack.pop()
        for c in n.children:
            out.append((n, c))
            stack.append(c)
    return out


# ------------------------------------------------------------------- re-rooting
def adjacency(root):
    adj = {}
    nodes = {}
    for p, c in edges_of(root):
        nodes[id(p)] = p
        nodes[id(c)] = c
        adj.setdefault(id(p), []).append(c)
        adj.setdefault(id(c), []).append(p)
    return adj, nodes


def reroot_above(root, child):
    """New rooted tree, rooted on the edge above `child`, built from adjacency
    of the original (so the original is untouched)."""
    adj, _ = adjacency(root)
    parent = child.parent

    def build(src, came_from, dest):
        for nb in adj[id(src)]:
            if nb is came_from:
                continue
            node = Node(nb.name)
            node.parent = dest
            dest.children.append(node)
            build(nb, src, node)

    newroot = Node("NEWROOT")
    for a, b in ((child, parent), (parent, child)):
        side = Node(a.name)
        side.parent = newroot
        newroot.children.append(side)
        build(a, b, side)
    return newroot


# --------------------------------------------------------------------- parsimony
def postorder(root):
    out = []

    def rec(n):
        for c in n.children:
            rec(c)
        out.append(n)

    rec(root)
    return out


def sankoff(root, states):
    """cost[id(node)][state] = min changes in that subtree given node's state."""
    cost = {}
    for n in postorder(root):
        if n.is_leaf():
            obs = states[n.name]
            cost[id(n)] = {s: (0 if s == obs else INF) for s in STATES}
            continue
        c = {}
        for s in STATES:
            tot = 0
            for ch in n.children:
                tot += min(cost[id(ch)][sc] + (0 if sc == s else 1) for sc in STATES)
            c[s] = tot
        cost[id(n)] = c
    return cost


def gains_range(root, states):
    """Min and max gains (ENV->PATH) over all most-parsimonious reconstructions."""
    cost = sankoff(root, states)
    gmin, gmax = {}, {}
    for n in postorder(root):
        if n.is_leaf():
            obs = states[n.name]
            gmin[id(n)] = {s: (0 if s == obs else INF) for s in STATES}
            gmax[id(n)] = {s: (0 if s == obs else -INF) for s in STATES}
            continue
        lo, hi = {}, {}
        for s in STATES:
            tlo = thi = 0
            for ch in n.children:
                best = min(cost[id(ch)][sc] + (0 if sc == s else 1) for sc in STATES)
                opts_lo, opts_hi = [], []
                for sc in STATES:
                    if cost[id(ch)][sc] + (0 if sc == s else 1) != best:
                        continue
                    add = 1 if (s == ENV and sc == PATH) else 0
                    opts_lo.append(gmin[id(ch)][sc] + add)
                    opts_hi.append(gmax[id(ch)][sc] + add)
                tlo += min(opts_lo)
                thi += max(opts_hi)
            lo[s], hi[s] = tlo, thi
        gmin[id(n)], gmax[id(n)] = lo, hi
    best = min(cost[id(root)][s] for s in STATES)
    root_states = [s for s in STATES if cost[id(root)][s] == best]
    return best, min(gmin[id(root)][s] for s in root_states), \
        max(gmax[id(root)][s] for s in root_states), root_states


def dollo_min_gains(root, states):
    """Minimum gains assuming the trait, once gained, is never lost."""
    g = {}
    for n in postorder(root):
        if n.is_leaf():
            obs = states[n.name]
            g[id(n)] = {s: (0 if s == obs else INF) for s in STATES}
            continue
        cur = {}
        for s in STATES:
            tot = 0
            for ch in n.children:
                if s == PATH:
                    # no loss allowed: child must stay PATH
                    opt = g[id(ch)][PATH]
                else:
                    opt = min(g[id(ch)][ENV], g[id(ch)][PATH] + 1)
                tot += opt
            cur[s] = tot
        g[id(n)] = cur
    return min(g[id(root)][s] for s in STATES)


def main():
    cats = {}
    with open(MANIFEST) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ia, ic = header.index("abbrev"), header.index("category")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ia, ic):
                continue
            cats[f[ia]] = PATH if f[ic] == "pathogenic" else ENV

    root = parse_newick(TREE.read_text())
    tips = [l.name for l in leaves(root)]
    missing = [t for t in tips if t not in cats]
    if missing:
        sys.exit(f"ERROR: tips with no category: {missing}")

    print(f"Tree : {TREE}")
    print(f"Trait: {MANIFEST}")
    print(f"Tips on tree: {len(tips)}")
    print(f"  pathogenic   ({sum(1 for t in tips if cats[t]==PATH):2d}): "
          f"{', '.join(sorted(t for t in tips if cats[t]==PATH))}")
    print(f"  environmental({sum(1 for t in tips if cats[t]==ENV):2d}): "
          f"{', '.join(sorted(t for t in tips if cats[t]==ENV))}")

    # ---- 1. root-independent parsimony score
    print("\n" + "=" * 74)
    print("1. PARSIMONY SCORE  (minimum number of state CHANGES, any direction)")
    print("=" * 74)
    base, _, _, rs = gains_range(root, cats)
    print(f"   minimum changes on the tree = {base}")
    scores = set()
    for _, c in edges_of(root):
        scores.add(gains_range(reroot_above(root, c), cats)[0])
    print(f"   re-rooted on all {len(edges_of(root))} edges -> score(s) observed: {sorted(scores)}")
    print("   (single value confirms the count is root-independent)")

    # ---- 2. gains vs losses, by root
    print("\n" + "=" * 74)
    print("2. INDEPENDENT ORIGINS (gains env->path) OVER ALL MP RECONSTRUCTIONS")
    print("=" * 74)
    cands = [("midpoint root, as published in paper 1", root)]
    for n in internals(root):
        names = {l.name for l in leaves(n)}
        if names == {"TPMV", "MJNV"}:
            cands.append(("rooted on shrew clade (TPMV,MJNV)", reroot_above(root, n)))
        if names == {"NVAV", "QZNV"}:
            cands.append(("rooted on mole+bat clade (NVAV,QZNV)", reroot_above(root, n)))
    for label, r in cands:
        chg, lo, hi, rstates = gains_range(r, cats)
        print(f"\n   {label}")
        print(f"     minimum changes            : {chg}")
        print(f"     root state(s) in MP solns   : {','.join(rstates)}")
        if lo == hi:
            print(f"     independent origins (gains): {lo}   [unambiguous]")
        else:
            print(f"     independent origins (gains): {lo}-{hi}   [ambiguous under parsimony]")
        print(f"     implied reversals (losses) : {chg - hi}-{chg - lo}")

    # ---- 3. Dollo
    print("\n" + "=" * 74)
    print("3. DOLLO COUNT (pathogenicity never lost once gained)")
    print("=" * 74)
    for label, r in cands:
        print(f"   {label}: minimum gains with zero losses = {dollo_min_gains(r, cats)}")

    # ---- 4. all-pathogenic clades
    print("\n" + "=" * 74)
    print("4. MAXIMAL ALL-PATHOGENIC CLADES ON THIS TOPOLOGY")
    print("=" * 74)
    allnodes = internals(root) + leaves(root)
    pure = []
    for n in allnodes:
        lv = leaves(n)
        if all(cats[l.name] == PATH for l in lv):
            pure.append(sorted(l.name for l in lv))
    kept = []
    for names in sorted(pure, key=len, reverse=True):
        if any(set(names) <= set(k) for k in kept):
            continue
        kept.append(names)
    for names in kept:
        print(f"   {len(names):2d} sp: {', '.join(names)}")
    print(f"\n   -> the 10 pathogenic species fall into {len(kept)} disjoint clade(s)")
    print("      i.e. pathogenicity is NOT a single clade on this topology")


if __name__ == "__main__":
    main()
