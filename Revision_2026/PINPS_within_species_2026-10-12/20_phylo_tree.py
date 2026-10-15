"""Minimal newick / rooted-tree / VCV utilities. Hand-rolled: no ete3, no dendropy.

Used by 21_phylo_stats.py. Pure numpy. Deterministic.
"""
import math

import numpy as np


# ---------------------------------------------------------------- newick parse
class Node:
    __slots__ = ("name", "blen", "children", "parent")

    def __init__(self, name="", blen=0.0):
        self.name = name
        self.blen = blen
        self.children = []
        self.parent = None

    def is_tip(self):
        return not self.children


def parse_newick(text):
    """Parse a newick string into a Node tree. Top node may be multifurcating."""
    s = text.strip()
    s = s[: s.rindex(";")] if ";" in s else s
    i = 0

    def parse_node():
        nonlocal i
        node = Node()
        if s[i] == "(":
            i += 1
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                if s[i] == ",":
                    i += 1
                    continue
                if s[i] == ")":
                    i += 1
                    break
                raise ValueError("bad newick at %d: %r" % (i, s[i : i + 20]))
        # label
        j = i
        while i < len(s) and s[i] not in "(),:;":
            i += 1
        node.name = s[j:i].strip()
        if i < len(s) and s[i] == ":":
            i += 1
            j = i
            while i < len(s) and s[i] not in "(),;":
                i += 1
            node.blen = float(s[j:i])
        return node

    root = parse_node()
    return root


def tips(node):
    out = []
    stack = [node]
    while stack:
        n = stack.pop()
        if n.is_tip():
            out.append(n)
        else:
            stack.extend(n.children)
    return out


def tip_names(node):
    return sorted(t.name for t in tips(node))


# ---------------------------------------------------------------- pruning
def prune(node, keep):
    """Return a copy of the subtree containing only tips in `keep`,
    collapsing degree-2 internal nodes and summing their branch lengths."""
    if node.is_tip():
        if node.name in keep:
            n = Node(node.name, node.blen)
            return n
        return None
    kids = [prune(c, keep) for c in node.children]
    kids = [k for k in kids if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        k = kids[0]
        k.blen += node.blen
        return k
    n = Node(node.name, node.blen)
    for k in kids:
        k.parent = n
        n.children.append(k)
    return n


# ---------------------------------------------------------------- adjacency
def adjacency(root):
    """Undirected adjacency: {node_id: [(other_id, length), ...]}, plus id->name."""
    nodes = []
    stack = [root]
    while stack:
        n = stack.pop()
        nodes.append(n)
        stack.extend(n.children)
    idx = {id(n): k for k, n in enumerate(nodes)}
    adj = {k: [] for k in range(len(nodes))}
    for n in nodes:
        for c in n.children:
            a, b = idx[id(n)], idx[id(c)]
            adj[a].append((b, c.blen))
            adj[b].append((a, c.blen))
    names = {idx[id(n)]: n.name for n in nodes}
    return adj, names


def _dists_from(adj, start):
    d = {start: 0.0}
    stack = [start]
    prev = {start: None}
    while stack:
        u = stack.pop()
        for v, L in adj[u]:
            if v not in d:
                d[v] = d[u] + L
                prev[v] = u
                stack.append(v)
    return d, prev


def rooted_from_point(adj, names, tip_order, edge, frac_from_u):
    """Build depth / root-to-tip path structures for a root placed on `edge`
    = (u, v, L) at distance frac_from_u*L from u.  If edge is None, the root is
    an existing node id given by frac_from_u (an int).

    Returns (depth_to_tip, paths) where paths[tip] is the ordered list of
    cumulative depths of ancestral nodes of that tip, root first.
    """
    # Build an augmented adjacency with a pseudo-root node if needed.
    adj2 = {k: list(v) for k, v in adj.items()}
    if edge is None:
        root = int(frac_from_u)
    else:
        u, v, L = edge
        root = max(adj2) + 1
        adj2[root] = []
        adj2[u] = [(w, l) for (w, l) in adj2[u] if w != v]
        adj2[v] = [(w, l) for (w, l) in adj2[v] if w != u]
        du = frac_from_u * L
        dv = L - du
        adj2[root].append((u, du))
        adj2[u].append((root, du))
        adj2[root].append((v, dv))
        adj2[v].append((root, dv))

    # DFS from root, recording node path and depth
    depth = {root: 0.0}
    parent = {root: None}
    order = [root]
    stack = [root]
    while stack:
        u = stack.pop()
        for v, L in adj2[u]:
            if v not in depth:
                depth[v] = depth[u] + L
                parent[v] = u
                order.append(v)
                stack.append(v)

    name_to_id = {}
    for k, nm in names.items():
        if nm:
            name_to_id[nm] = k
    paths = {}
    for t in tip_order:
        k = name_to_id[t]
        chain = []
        while k is not None:
            chain.append(k)
            k = parent[k]
        chain.reverse()
        paths[t] = chain
    return depth, paths


def vcv_rooted_at_edge(adj, names, tip_order, u, v, L, frac):
    depth, paths = rooted_from_point(adj, names, tip_order, (u, v, L), frac)
    n = len(tip_order)
    C = np.zeros((n, n))
    for i, a in enumerate(tip_order):
        for j, b in enumerate(tip_order):
            pa, pb = paths[a], paths[b]
            m = 0
            while m < len(pa) and m < len(pb) and pa[m] == pb[m]:
                m += 1
            C[i, j] = depth[pa[m - 1]]
    return C


def vcv_rooted_at_node(adj, names, tip_order, node_id):
    depth, paths = rooted_from_point(adj, names, tip_order, None, node_id)
    n = len(tip_order)
    C = np.zeros((n, n))
    for i, a in enumerate(tip_order):
        for j, b in enumerate(tip_order):
            pa, pb = paths[a], paths[b]
            m = 0
            while m < len(pa) and m < len(pb) and pa[m] == pb[m]:
                m += 1
            C[i, j] = depth[pa[m - 1]]
    return C


def midpoint_edge(adj, names, tip_order):
    """Return (u, v, L, frac) for the midpoint root of the tree."""
    name_to_id = {nm: k for k, nm in names.items() if nm}
    tid = [name_to_id[t] for t in tip_order]
    best = (-1.0, None, None)
    for a in tid:
        d, _ = _dists_from(adj, a)
        for b in tid:
            if d.get(b, -1) > best[0]:
                best = (d[b], a, b)
    total, a, b = best
    d, prev = _dists_from(adj, a)
    # walk from b back to a collecting the path
    chain = [b]
    while prev[chain[-1]] is not None:
        chain.append(prev[chain[-1]])
    chain.reverse()  # a ... b
    half = total / 2.0
    acc = 0.0
    for k in range(len(chain) - 1):
        u, v = chain[k], chain[k + 1]
        L = next(l for (w, l) in adj[u] if w == v)
        if acc + L >= half:
            frac = (half - acc) / L
            return u, v, L, frac, total
        acc += L
    raise RuntimeError("midpoint not found")


def find_edge_subtending(root, taxa):
    """In a rooted Node tree, find the node whose tip set == taxa (a clade), and
    return (parent_name_free) ids for use with adjacency: we return the Node."""
    want = set(taxa)
    stack = [root]
    while stack:
        n = stack.pop()
        if set(tip_names(n)) == want:
            return n
        stack.extend(n.children)
    return None


def to_newick(node):
    if node.is_tip():
        return "%s:%.10f" % (node.name, node.blen)
    inner = ",".join(to_newick(c) for c in node.children)
    return "(%s):%.10f" % (inner, node.blen)


# ------------------------------------------------- parsimony (Sankoff DP, exact)
def parsimony_score(root, states, nstates=2):
    """Minimum number of state changes for a discrete character on the tree.

    Exact Sankoff dynamic programming with unit change cost, so it is correct for
    multifurcating nodes too (plain Fitch is not).  With symmetric unit costs the
    score is independent of root placement.  `states` maps tip name -> int state.
    """

    def rec(n):
        if n.is_tip():
            return [0.0 if s == states[n.name] else math.inf for s in range(nstates)]
        tot = [0.0] * nstates
        for c in n.children:
            cc = rec(c)
            for s in range(nstates):
                tot[s] += min(cc[t] + (0 if t == s else 1) for t in range(nstates))
        return tot

    return int(min(rec(root)))
