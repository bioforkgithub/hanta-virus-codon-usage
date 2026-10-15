#!/usr/bin/env python3
"""Confound-direct statistical treatment of the per-species within-species piN/piS.

Attacks the four confounds head on at the species level (k = 10 valid species):
  (a) study effort   - regress piN/piS on available-sequence counts and on species
                       age (years since first GenBank deposit), then test the
                       disease effect on the residuals
  (b) demography     - test piS alone between groups
  (c) host ecology   - rodent vs non-rodent, tested directly
  (d) clade          - New World HPS clade membership (defined from the ML tree,
                       not from memory) as a predictor
and ranks every single predictor by variance explained, with and without a
phylogenetic control against S22_ml_clean.nwk.

Everything is hand-rolled on numpy (exact permutation nulls where k allows),
cross-checked against scipy where scipy provides the same quantity.
Deterministic: master seed 20261012, no wall-clock input.

Usage:
  python3 -I 21_stats_confound_direct.py <ANALYSIS_ROOT> <LABELS_TSV> <TREE_NWK> <OUTDIR>
"""
import csv
import itertools
import json
import math
import os
import sys

import numpy as np

SEED = 20261012
TWO_SIDED_ALPHA = 0.05

# --------------------------------------------------------------------------
# basic hand-rolled statistics
# --------------------------------------------------------------------------


def rankdata(a):
    """Mid-ranks (average of tied ranks), 1-based."""
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)
    sa = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and sa[j + 1] == sa[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


_MWU_NULL = {}


def _mwu_null(r, k1):
    """Exact permutation distribution of U1 for the given mid-rank vector.

    When there are no ties the rank vector is always 1..k, so the null depends
    only on (k, k1) and is cached; with ties the full enumeration is redone.
    """
    k = len(r)
    tied = len(set(r)) < k
    key = (k, k1, tuple(np.round(r, 9)) if tied else None)
    if key in _MWU_NULL:
        return _MWU_NULL[key]
    off = k1 * (k1 + 1) / 2.0
    vals = np.fromiter((sum(r[i] for i in comb) - off
                        for comb in itertools.combinations(range(k), k1)),
                       dtype=float, count=math.comb(k, k1))
    _MWU_NULL[key] = vals
    return vals


def mwu_exact(y, g):
    """Exact two-sided Mann-Whitney U by full enumeration of group assignments.

    g is a boolean array; group 1 = g True.  Mid-ranks are used, and the null is
    the exact permutation distribution of the group-1 rank sum, so ties are
    handled exactly.  Returns a dict.
    """
    y = np.asarray(y, dtype=float)
    g = np.asarray(g, dtype=bool)
    k = len(y)
    k1 = int(g.sum())
    k2 = k - k1
    if k1 == 0 or k2 == 0:
        return {"n1": k1, "n2": k2, "U1": None, "p_two_sided": None}
    r = rankdata(y)
    R1 = r[g].sum()
    U1 = R1 - k1 * (k1 + 1) / 2.0
    EU = k1 * k2 / 2.0
    dev = abs(U1 - EU)
    null = _mwu_null(r, k1)
    n_tot = len(null)
    n_ge = int(np.count_nonzero(np.abs(null - EU) >= dev - 1e-12))
    # Hodges-Lehmann shift (median of all cross-group differences, g=True minus g=False)
    diffs = np.subtract.outer(y[g], y[~g]).ravel()
    return {"n1": k1, "n2": k2,
            "U1": float(U1), "EU": float(EU),
            "rank_sum_group1": float(R1),
            "rank_biserial_r": float(2.0 * U1 / (k1 * k2) - 1.0),
            "p_two_sided": float(n_ge / n_tot),
            "n_assignments": int(n_tot),
            "p_floor_two_sided": float(2.0 / n_tot),
            "median_group1": float(np.median(y[g])),
            "median_group2": float(np.median(y[~g])),
            "mean_group1": float(np.mean(y[g])),
            "mean_group2": float(np.mean(y[~g])),
            "hodges_lehmann_g1_minus_g2": float(np.median(diffs))}


_SPEARMAN_NULL = {}


def spearman_null_exact(k):
    """Exact null distribution of sum(i * perm(i)) for ranks 1..k, no ties.

    With no ties every variable has rank vector 1..k, so one cached null serves
    every variable pair at that k.  Enumerated exactly for k <= 10.
    """
    if k in _SPEARMAN_NULL:
        return _SPEARMAN_NULL[k]
    base = np.arange(1, k + 1, dtype=np.int64)
    n = math.factorial(k)
    out = np.empty(n, dtype=np.int64)
    it = itertools.permutations(range(1, k + 1))
    pos = 0
    chunk = 200000
    while pos < n:
        take = min(chunk, n - pos)
        blk = np.fromiter(itertools.chain.from_iterable(
            itertools.islice(it, take)), dtype=np.int64,
            count=take * k).reshape(take, k)
        out[pos:pos + take] = blk @ base
        pos += take
    out.sort()
    vals = out.astype(float)
    _SPEARMAN_NULL[k] = vals
    return vals


def _perm_matrix(k, n, rng):
    """n independent random permutations of 0..k-1, one per row."""
    return np.argsort(rng.random((n, k)), axis=1)


def spearman(x, y, n_mc=19999, seed=SEED):
    """Spearman rho with an exact permutation p where possible, else Monte Carlo."""
    rx, ry = rankdata(x), rankdata(y)
    k = len(rx)
    sx, sy = rx.std(), ry.std()
    if sx == 0 or sy == 0:
        return {"rho": None, "p_two_sided": None, "method": "degenerate", "k": k}
    rho = float(np.corrcoef(rx, ry)[0, 1])
    ties = (len(set(rx)) < k) or (len(set(ry)) < k)
    if k <= 10 and not ties:
        S = float(np.dot(rx, ry))
        null = spearman_null_exact(k)
        mu = null.mean()
        dev = abs(S - mu)
        p = float(np.mean(np.abs(null - mu) >= dev - 1e-9))
        method = "exact_enumeration_%d!" % k
        n_null = len(null)
    else:
        rng = np.random.default_rng([seed, 7, k])
        P = _perm_matrix(k, n_mc, rng)
        a = rx - rx.mean()
        b = ry - ry.mean()
        num = (ry[P] - ry.mean()) @ a
        den = math.sqrt(float(a @ a)) * math.sqrt(float(b @ b))
        null = num / den
        p = float((np.count_nonzero(np.abs(null) >= abs(rho) - 1e-12) + 1)
                  / (n_mc + 1))
        method = "monte_carlo_%d" % n_mc
        n_null = n_mc
    return {"rho": rho, "rho_sq": rho * rho, "p_two_sided": p,
            "method": method, "n_null": n_null, "k": k}


def ols(y, X, names):
    """OLS with an intercept already in X.  Returns coefficients, t, R2, AICc."""
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    k, p = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    rss = float(resid @ resid)
    tss = float(((y - y.mean()) ** 2).sum())
    dof = k - p
    out = {"k": k, "n_params": p, "rss": rss, "tss": tss,
           "r2": (1.0 - rss / tss) if tss > 0 else None,
           "r2_adj": (1.0 - (rss / dof) / (tss / (k - 1)))
           if (tss > 0 and dof > 0) else None,
           "resid": resid.tolist(), "coef": {}, "dof": dof}
    if dof > 0 and rss > 0:
        s2 = rss / dof
        XtXi = np.linalg.pinv(X.T @ X)
        se = np.sqrt(np.diag(XtXi) * s2)
        for i, nm in enumerate(names):
            t = beta[i] / se[i] if se[i] > 0 else None
            out["coef"][nm] = {"beta": float(beta[i]), "se": float(se[i]),
                               "t": (float(t) if t is not None else None),
                               "p_two_sided": (student_p(t, dof)
                                               if t is not None else None)}
        ll = -0.5 * k * (math.log(2 * math.pi * rss / k) + 1.0)
        npar = p + 1
        aic = -2 * ll + 2 * npar
        out["loglik"] = float(ll)
        out["aic"] = float(aic)
        out["aicc"] = float(aic + (2 * npar * (npar + 1)) / (k - npar - 1)) \
            if k - npar - 1 > 0 else None
    else:
        for i, nm in enumerate(names):
            out["coef"][nm] = {"beta": float(beta[i]), "se": None,
                               "t": None, "p_two_sided": None}
    return out


def student_p(t, dof):
    """Two-sided Student-t p-value by continued-fraction incomplete beta."""
    t = abs(float(t))
    x = dof / (dof + t * t)
    return float(betainc(0.5 * dof, 0.5, x))


def betainc(a, b, x):
    """Regularised incomplete beta I_x(a,b) (Lentz continued fraction)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1) / (a + b + 2):
        return math.exp(lbeta) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbeta) * _betacf(b, a, 1 - x) / b


def _betacf(a, b, x, itmax=300, eps=3e-16):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < 1e-300:
        d = 1e-300
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h


def kruskal(groups):
    """Kruskal-Wallis H with tie correction (same convention as the R2.2 scripts)."""
    allv = np.concatenate(groups)
    k = len(allv)
    r = rankdata(allv)
    pos = 0
    H = 0.0
    for gset in groups:
        n = len(gset)
        Ri = r[pos:pos + n].sum()
        H += Ri * Ri / n
        pos += n
    H = 12.0 / (k * (k + 1)) * H - 3.0 * (k + 1)
    _, cnt = np.unique(allv, return_counts=True)
    tie = 1.0 - sum(c ** 3 - c for c in cnt) / float(k ** 3 - k)
    if tie > 0:
        H /= tie
    df = len(groups) - 1
    return {"H": float(H), "df": df,
            "eta_sq_rank": float((H - df) / (k - df - 1)) if k - df - 1 > 0 else None,
            "n_total": k, "group_sizes": [len(g) for g in groups]}


def kruskal_perm_p(groups, n_perm=99999, seed=SEED):
    """Exact (enumerated) or Monte-Carlo permutation p for Kruskal-Wallis H."""
    sizes = [len(g) for g in groups]
    allv = np.concatenate(groups)
    k = len(allv)
    obs = kruskal(groups)["H"]
    n_assign = 1
    rem = k
    for s in sizes[:-1]:
        n_assign *= math.comb(rem, s)
        rem -= s
    if n_assign <= 200000:
        cnt = tot = 0
        for assign in _assignments(list(range(k)), sizes):
            gs = [allv[list(a)] for a in assign]
            if kruskal(gs)["H"] >= obs - 1e-12:
                cnt += 1
            tot += 1
        return {"p_two_sided": float(cnt / tot), "method": "exact_%d" % tot,
                "n_null": tot}
    rng = np.random.default_rng([seed, 11, k])
    cnt = 0
    for _ in range(n_perm):
        pr = rng.permutation(allv)
        gs, pos = [], 0
        for s in sizes:
            gs.append(pr[pos:pos + s])
            pos += s
        if kruskal(gs)["H"] >= obs - 1e-12:
            cnt += 1
    return {"p_two_sided": float((cnt + 1) / (n_perm + 1)),
            "method": "monte_carlo_%d" % n_perm, "n_null": n_perm}


def _assignments(idx, sizes):
    if len(sizes) == 1:
        yield [tuple(idx)]
        return
    for comb in itertools.combinations(idx, sizes[0]):
        rest = [i for i in idx if i not in comb]
        for tail in _assignments(rest, sizes[1:]):
            yield [comb] + tail


# --------------------------------------------------------------------------
# phylogenetics
# --------------------------------------------------------------------------


def read_tree(path):
    from Bio import Phylo
    return Phylo.read(path, "newick")


def patristic(tree, tips):
    k = len(tips)
    D = np.zeros((k, k))
    nodes = {t.name: t for t in tree.get_terminals()}
    for i in range(k):
        for j in range(i + 1, k):
            d = tree.distance(nodes[tips[i]], nodes[tips[j]])
            D[i, j] = D[j, i] = d
    return D


def vcv_from_rooted(tree, tips):
    """BM variance-covariance: C_ij = root-to-MRCA(i,j) path length."""
    k = len(tips)
    nodes = {t.name: t for t in tree.get_terminals()}
    root = tree.root
    C = np.zeros((k, k))
    for i in range(k):
        C[i, i] = tree.distance(root, nodes[tips[i]])
        for j in range(i + 1, k):
            m = tree.common_ancestor([nodes[tips[i]], nodes[tips[j]]])
            C[i, j] = C[j, i] = tree.distance(root, m)
    return C


def pgls(y, X, C, names):
    """GLS under a BM covariance C, with t-tests and a GLS R-squared."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    k, p = X.shape
    L = np.linalg.cholesky(C + 1e-12 * np.eye(k))
    Li = np.linalg.inv(L)
    yw = Li @ y
    Xw = Li @ X
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    res = yw - Xw @ beta
    rss = float(res @ res)
    one = Li @ np.ones((k, 1))
    b0, *_ = np.linalg.lstsq(one, yw, rcond=None)
    tss = float(((yw - one.ravel() * b0[0]) ** 2).sum())
    dof = k - p
    out = {"k": k, "rss_gls": rss, "tss_gls": tss,
           "r2_gls": (1.0 - rss / tss) if tss > 0 else None,
           "dof": dof, "coef": {}}
    if dof > 0 and rss > 0:
        s2 = rss / dof
        XtXi = np.linalg.pinv(Xw.T @ Xw)
        se = np.sqrt(np.diag(XtXi) * s2)
        for i, nm in enumerate(names):
            t = beta[i] / se[i] if se[i] > 0 else None
            out["coef"][nm] = {"beta": float(beta[i]), "se": float(se[i]),
                               "t": float(t) if t is not None else None,
                               "p_two_sided": student_p(t, dof)
                               if t is not None else None}
        ll = -0.5 * (k * math.log(2 * math.pi * rss / k) + k) \
            - float(np.log(np.diag(L)).sum())
        npar = p + 1
        aic = -2 * ll + 2 * npar
        out["loglik"] = float(ll)
        out["aic"] = float(aic)
        out["aicc"] = float(aic + (2 * npar * (npar + 1)) / (k - npar - 1)) \
            if k - npar - 1 > 0 else None
    return out


def _blomberg_rows(Y, Ci, one, expect):
    """Blomberg's K for each ROW of Y (vectorised over permutations)."""
    a = (Y @ (Ci @ one)) / float(one @ Ci @ one)
    D = Y - a[:, None]
    mse0 = np.einsum("ij,ij->i", D, D)
    mse = np.einsum("ij,jk,ik->i", D, Ci, D)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (mse0 / mse) / expect


def blomberg_k(y, C, n_perm=9999, seed=SEED):
    """Blomberg et al. (2003) K, with a tip-label permutation p (K>1 = clumped)."""
    y = np.asarray(y, float)
    k = C.shape[0]
    Ci = np.linalg.inv(C)
    one = np.ones(k)
    expect = (np.trace(C) - k / float(one @ Ci @ one)) / (k - 1)
    obs = float(_blomberg_rows(y[None, :], Ci, one, expect)[0])
    rng = np.random.default_rng([seed, 23, k])
    null = _blomberg_rows(y[_perm_matrix(k, n_perm, rng)], Ci, one, expect)
    p = float((np.count_nonzero(null >= obs - 1e-12) + 1) / (n_perm + 1))
    return {"blomberg_K": obs, "perm_p_greater": p, "n_perm": n_perm,
            "null_median": float(np.median(null))}


def mantel(D, y, n_perm=9999, seed=SEED):
    """Mantel test: patristic distance vs |trait difference|.  Rooting-invariant."""
    k = len(y)
    iu = np.triu_indices(k, 1)
    d = D[iu]
    dz = d - d.mean()
    dn = math.sqrt(float(dz @ dz))
    y = np.asarray(y, float)

    def rvals(Y):
        T = np.abs(Y[:, iu[0]] - Y[:, iu[1]])
        Tz = T - T.mean(axis=1, keepdims=True)
        nn = np.sqrt(np.einsum("ij,ij->i", Tz, Tz))
        with np.errstate(divide="ignore", invalid="ignore"):
            return (Tz @ dz) / (nn * dn)
    obs = float(rvals(y[None, :])[0])
    rng = np.random.default_rng([seed, 31, k])
    null = rvals(y[_perm_matrix(k, n_perm, rng)])
    p = float((np.count_nonzero(np.abs(null) >= abs(obs) - 1e-12) + 1)
              / (n_perm + 1))
    return {"mantel_r": obs, "p_two_sided": p, "n_perm": n_perm}


def morans_i(D, y, n_perm=9999, seed=SEED):
    """Moran's I with w_ij = 1/patristic distance.  Rooting-invariant."""
    k = len(y)
    W = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            if i != j and D[i, j] > 0:
                W[i, j] = 1.0 / D[i, j]
    S0 = W.sum()
    y = np.asarray(y, float)

    def ivals(Y):
        Z = Y - Y.mean(axis=1, keepdims=True)
        den = np.einsum("ij,ij->i", Z, Z)
        num = np.einsum("ij,jk,ik->i", Z, W, Z)
        with np.errstate(divide="ignore", invalid="ignore"):
            return (k / S0) * num / den
    obs = float(ivals(y[None, :])[0])
    rng = np.random.default_rng([seed, 37, k])
    null = ivals(y[_perm_matrix(k, n_perm, rng)])
    p = float((np.count_nonzero(np.abs(null) >= abs(obs) - 1e-12) + 1)
              / (n_perm + 1))
    return {"morans_I": obs, "expected_I": -1.0 / (k - 1),
            "p_two_sided": p, "n_perm": n_perm}


def phylo_anova(y, g, C, n_sim=10000, seed=SEED):
    """Garland et al. (1993) phylogenetic ANOVA: observed F vs BM-simulated F."""
    y = np.asarray(y, float)
    g = np.asarray(g, bool)
    k = len(y)

    n1, n2 = int(g.sum()), k - int(g.sum())

    def fstat_rows(Y):
        """One-way ANOVA F for each row of Y."""
        m1 = Y[:, g].mean(axis=1)
        m2 = Y[:, ~g].mean(axis=1)
        gm = Y.mean(axis=1)
        ssb = n1 * (m1 - gm) ** 2 + n2 * (m2 - gm) ** 2
        ssw = (np.einsum("ij,ij->i", Y[:, g] - m1[:, None],
                         Y[:, g] - m1[:, None])
               + np.einsum("ij,ij->i", Y[:, ~g] - m2[:, None],
                           Y[:, ~g] - m2[:, None]))
        with np.errstate(divide="ignore", invalid="ignore"):
            f = ssb / (ssw / (k - 2))
        return np.where(ssw <= 0, np.inf, f)
    obs = float(fstat_rows(y[None, :])[0])
    # BM parameters by ML under the tree
    Ci = np.linalg.inv(C)
    one = np.ones(k)
    a = float((one @ Ci @ y) / (one @ Ci @ one))
    d = y - a
    sig2 = float(d @ Ci @ d) / k
    L = np.linalg.cholesky(C + 1e-12 * np.eye(k))
    rng = np.random.default_rng([seed, 41, k])
    Z = rng.standard_normal((n_sim, k))
    sims = fstat_rows(a + math.sqrt(sig2) * (Z @ L.T))
    cnt = int(np.count_nonzero(sims >= obs - 1e-12))
    # the ordinary (star-tree) p for comparison
    p_star = student_p(math.sqrt(obs) if obs != np.inf else 1e9, k - 2) \
        if obs != np.inf else 0.0
    return {"F_obs": obs, "p_phylo_anova": float((cnt + 1) / (n_sim + 1)),
            "p_ordinary_anova": p_star, "n_sim": n_sim,
            "bm_root_state": a, "bm_sigma2": sig2,
            "F_null_95th": float(np.percentile(sims, 95))}


# --------------------------------------------------------------------------
# data assembly
# --------------------------------------------------------------------------


def load_labels(path):
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            out[f["abbrev"]] = f
    return out


def load_table(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(s):
    s = (s or "").strip()
    if s in ("", "NA", "nan", "None"):
        return None
    return float(s)


def build_design(root, labels, tree, effort, table):
    """Per-species predictor table (one row per species, set-independent parts)."""
    nw_clade = sorted(x.name for x in
                      tree.common_ancestor(["SNV", "BAYV"]).get_terminals())
    whole = {r["abbrev"]: r for r in table
             if r["seqset"] == "ss_all" and r["n_level"] == "full"}
    n_usable = {}
    for ab in labels:
        p = os.path.join(root, "logs", ab + ".summary.json")
        if os.path.exists(p):
            with open(p) as fh:
                n_usable[ab] = json.load(fh).get("n_usable")
    meta_path = os.path.join(root, "results", "sets", "sequence_metadata.tsv")
    first_year = {}
    with open(meta_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            y = r["year"].strip()
            if y:
                ab = r["abbrev"]
                yi = int(y)
                if ab not in first_year or yi < first_year[ab]:
                    first_year[ab] = yi
    design = {}
    for ab, lab in labels.items():
        es = effort["species"].get(ab, {})
        design[ab] = {
            "abbrev": ab,
            "species": lab["species"],
            "category": lab["category"],
            "category_tulv_moved": ("pathogenic" if ab == "TULV"
                                    else lab["category"]),
            "host_class": lab["host_class"],
            "is_rodent": lab["is_rodent"] == "True",
            "nw_hps_clade": ab in nw_clade,
            "esearch_count": es.get("esearch_count_1975_2026"),
            "first_genbank_year": es.get("first_genbank_year"),
            "n_usable_records": n_usable.get(ab),
            "n_whole_species_ss": (int(whole[ab]["n_sequences"])
                                   if ab in whole else None),
            "first_collection_year": first_year.get(ab),
        }
    return design, nw_clade


BIN_PREDICTORS = [
    ("cat_pathogenic", "category", "pathogenic", "(a') disease association"),
    ("cat_pathogenic_TULVmoved", "category_tulv_moved", "pathogenic",
     "(a') disease, TULV relabelled pathogenic (R3.3c sensitivity)"),
    ("is_rodent", "is_rodent", True, "(c) host ecology: rodent reservoir"),
    ("nw_hps_clade", "nw_hps_clade", True,
     "(d) clade: New World HPS clade (MRCA(SNV,BAYV) of S22_ml_clean.nwk)"),
]

CONT_PREDICTORS = [
    ("log10_esearch_count", "(a) study effort: log10 GenBank records matching the "
     "download query (NCBI, all 22 species)"),
    ("log10_n_whole_species_ss", "(a) study effort: log10 usable complete N CDS, "
     "whole species sensu stricto"),
    ("years_since_first_genbank", "(a) species age proxy: 2026 minus the year of "
     "the first GenBank deposit"),
    ("piS", "(b) demography: synonymous diversity of the same set"),
    ("n_countries", "sampling breadth of the valid set"),
    ("frac_human", "fraction of the valid set that is human-derived"),
]


def dataset(table, design, seqset, n_level, require_usable, min_n=10):
    # n_countries / n_human are recorded only on the `full` row of each set, so
    # subsample rows inherit them from their own parent set.  They therefore
    # describe the PARENT single-population set, not the particular draw.
    parent = {r["abbrev"]: r for r in table
              if r["seqset"] == seqset and r["n_level"] == "full"}
    rows = []
    for r in table:
        if r["seqset"] != seqset or r["n_level"] != n_level:
            continue
        if require_usable and r["usable_for_group_test"] != "True":
            continue
        if int(r["n_sequences"]) < min_n:
            continue
        ab = r["abbrev"]
        d = dict(design[ab])
        d["n_sequences"] = int(r["n_sequences"])
        par = parent.get(ab, {})
        d["n_countries"] = fnum(r["n_countries"]) or fnum(par.get("n_countries"))
        d["n_human"] = (fnum(r["n_human"]) if fnum(r["n_human"]) is not None
                        else fnum(par.get("n_human")))
        n_par = int(par["n_sequences"]) if par else int(r["n_sequences"])
        d["frac_human"] = (d["n_human"] / n_par
                           if d["n_human"] is not None else None)
        d["n_parent_set"] = n_par
        d["piN_piS"] = fnum(r["piN_piS"])
        d["piS"] = fnum(r["piS"])
        d["piN"] = fnum(r["piN"])
        d["piN_piS_jc"] = fnum(r["piN_piS_jc"])
        d["flag_piS_saturated"] = r["flag_piS_saturated"] == "True"
        d["log10_esearch_count"] = (math.log10(d["esearch_count"])
                                    if d["esearch_count"] else None)
        d["log10_n_whole_species_ss"] = (math.log10(d["n_whole_species_ss"])
                                         if d["n_whole_species_ss"] else None)
        d["years_since_first_genbank"] = (2026 - d["first_genbank_year"]
                                          if d["first_genbank_year"] else None)
        rows.append(d)
    rows.sort(key=lambda d: d["abbrev"])
    return rows


# --------------------------------------------------------------------------
# the battery
# --------------------------------------------------------------------------

RESPONSES = ["piN_piS", "piS", "piN"]


def run_battery(rows, dsname, note, tree, out_tests, out_json):
    tips = [d["abbrev"] for d in rows]
    k = len(rows)
    rec = {"dataset": dsname, "note": note, "k": k, "species": tips,
           "n_sequences_used": {d["abbrev"]: d["n_sequences"] for d in rows},
           "saturated_flagged": [d["abbrev"] for d in rows
                                 if d["flag_piS_saturated"]]}

    def add(test, predictor, response, res, extra=None):
        row = {"dataset": dsname, "k": k, "test": test,
               "predictor": predictor, "response": response}
        row.update({kk: vv for kk, vv in res.items()
                    if kk not in ("resid",)})
        if extra:
            row.update(extra)
        out_tests.append(row)

    # ---- (1) binary predictors, exact Mann-Whitney, every response
    for pname, field, truth, desc in BIN_PREDICTORS:
        g = np.array([rows[i][field] == truth for i in range(k)], dtype=bool)
        if g.sum() == 0 or g.sum() == k:
            add("exact_MWU", pname, "ALL", {"skipped": "degenerate split",
                                            "n1": int(g.sum()), "n2": k - int(g.sum())})
            continue
        for resp in RESPONSES:
            y = np.array([d[resp] for d in rows], float)
            add("exact_MWU", pname, resp, mwu_exact(y, g),
                {"description": desc,
                 "group1_members": ",".join(sorted(t for t, gg in zip(tips, g) if gg)),
                 "group2_members": ",".join(sorted(t for t, gg in zip(tips, g) if not gg))})

    # ---- host class as a 3-level factor (Kruskal-Wallis)
    hc = sorted(set(d["host_class"] for d in rows))
    if len(hc) >= 2:
        for resp in RESPONSES:
            groups = [np.array([d[resp] for d in rows if d["host_class"] == h], float)
                      for h in hc]
            res = kruskal(groups)
            res.update(kruskal_perm_p(groups))
            add("kruskal_wallis", "host_class_%dlevel" % len(hc), resp, res,
                {"levels": ",".join("%s:%d" % (h, len(g)) for h, g in zip(hc, groups))})

    # ---- (2) continuous predictors: Spearman
    for pname, desc in CONT_PREDICTORS:
        x = [d.get(pname) for d in rows]
        if any(v is None for v in x):
            add("spearman", pname, "ALL",
                {"skipped": "missing values", "n_missing":
                 sum(1 for v in x if v is None)})
            continue
        for resp in RESPONSES:
            if pname == resp:
                continue
            y = [d[resp] for d in rows]
            add("spearman", pname, resp, spearman(np.array(x, float),
                                                  np.array(y, float)),
                {"description": desc})

    # ---- (3) variance explained: single-predictor OLS, raw and rank scale
    ve = []
    allpred = []
    for pname, field, truth, desc in BIN_PREDICTORS:
        g = np.array([rows[i][field] == truth for i in range(k)], float)
        if 0 < g.sum() < k:
            allpred.append((pname, g, desc))
    for pname, desc in CONT_PREDICTORS:
        x = [d.get(pname) for d in rows]
        if not any(v is None for v in x):
            allpred.append((pname, np.array(x, float), desc))
    for resp in RESPONSES:
        yraw = np.array([d[resp] for d in rows], float)
        yrank = rankdata(yraw)
        for pname, x, desc in allpred:
            if pname == resp:
                continue
            for scale, y in (("raw", yraw), ("rank", yrank)):
                X = np.column_stack([np.ones(k), x])
                r = ols(y, X, ["intercept", pname])
                ve.append({"dataset": dsname, "k": k, "response": resp,
                           "scale": scale, "predictor": pname,
                           "r2": r["r2"], "r2_adj": r["r2_adj"],
                           "beta": r["coef"][pname]["beta"],
                           "t": r["coef"][pname]["t"],
                           "p_two_sided": r["coef"][pname]["p_two_sided"],
                           "aicc": r.get("aicc"), "description": desc})
    rec["variance_explained"] = ve

    # ---- (4) residual tests: disease effect net of each confound
    resid_tests = []
    for resp in ("piN_piS",):
        yraw = np.array([d[resp] for d in rows], float)
        yrank = rankdata(yraw)
        gcat = np.array([d["category"] == "pathogenic" for d in rows], bool)
        gcat2 = np.array([d["category_tulv_moved"] == "pathogenic" for d in rows],
                         bool)
        for cname in ["log10_esearch_count", "log10_n_whole_species_ss",
                      "years_since_first_genbank", "piS", "n_countries",
                      "frac_human"]:
            x = [d.get(cname) for d in rows]
            if any(v is None for v in x):
                continue
            for scale, y in (("raw", yraw), ("rank", yrank)):
                X = np.column_stack([np.ones(k), np.array(x, float)])
                fit = ols(y, X, ["intercept", cname])
                r = np.array(fit["resid"], float)
                for gname, g in (("cat_pathogenic", gcat),
                                 ("cat_pathogenic_TULVmoved", gcat2)):
                    if not (0 < g.sum() < k):
                        continue
                    m = mwu_exact(r, g)
                    resid_tests.append({
                        "dataset": dsname, "k": k, "response": resp,
                        "scale": scale, "covariate_removed": cname,
                        "covariate_r2": fit["r2"],
                        "covariate_p": fit["coef"][cname]["p_two_sided"],
                        "group": gname, "U1": m["U1"], "n1": m["n1"],
                        "n2": m["n2"], "p_two_sided": m["p_two_sided"],
                        "p_floor": m["p_floor_two_sided"],
                        "rank_biserial_r": m["rank_biserial_r"],
                        "median_resid_path": m["median_group1"],
                        "median_resid_env": m["median_group2"]})
        # multiple covariates at once (effort + demography), reported with dof warning
        cols, nms = [np.ones(k)], ["intercept"]
        for cname in ["log10_esearch_count", "piS"]:
            cols.append(np.array([d[cname] for d in rows], float))
            nms.append(cname)
        cols.append(np.array([1.0 if d["category"] == "pathogenic" else 0.0
                              for d in rows]))
        nms.append("cat_pathogenic")
        fit = ols(yrank, np.column_stack(cols), nms)
        rec["multivariate_rank_model"] = {
            "formula": "rank(piN_piS) ~ log10_esearch_count + piS + cat_pathogenic",
            "k": k, "r2": fit["r2"], "r2_adj": fit["r2_adj"],
            "dof": fit["dof"], "coef": fit["coef"],
            "warning": "k=%d with 4 parameters: reported for completeness only" % k}
    rec["residual_tests"] = resid_tests

    # ---- (5) phylogenetic control
    import copy
    from Bio import Phylo
    phy = {}
    present = [t for t in tips]
    for rootmode in ("midpoint", "outgroup_HTNV_SEOV_DOBV"):
        t2 = copy.deepcopy(tree)
        try:
            if rootmode == "midpoint":
                t2.root_at_midpoint()
            else:
                t2.root_with_outgroup(["HTNV", "SEOV", "DOBV"])
        except Exception as exc:  # noqa: BLE001
            phy[rootmode] = {"error": str(exc)}
            continue
        for leaf in [x for x in t2.get_terminals() if x.name not in present]:
            t2.prune(leaf)
        tip_order = [x.name for x in t2.get_terminals()]
        C = vcv_from_rooted(t2, tip_order)
        idx = [tip_order.index(t) for t in present]
        C = C[np.ix_(idx, idx)]
        sub = {"root_mode": rootmode, "tip_order": present,
               "vcv_diag": np.diag(C).tolist()}
        for resp in RESPONSES:
            y = np.array([d[resp] for d in rows], float)
            entry = {}
            entry["phylogenetic_signal"] = blomberg_k(y, C)
            for pname, field, truth, desc in BIN_PREDICTORS:
                g = np.array([rows[i][field] == truth for i in range(k)], bool)
                if not (0 < g.sum() < k):
                    continue
                X = np.column_stack([np.ones(k), g.astype(float)])
                pg = pgls(y, X, C, ["intercept", pname])
                entry["pgls_" + pname] = {
                    "beta": pg["coef"][pname]["beta"],
                    "t": pg["coef"][pname]["t"],
                    "p_two_sided": pg["coef"][pname]["p_two_sided"],
                    "r2_gls": pg["r2_gls"], "aicc": pg.get("aicc")}
                entry["phylo_anova_" + pname] = phylo_anova(y, g, C)
            for pname, desc in CONT_PREDICTORS:
                x = [d.get(pname) for d in rows]
                if any(v is None for v in x) or pname == resp:
                    continue
                X = np.column_stack([np.ones(k), np.array(x, float)])
                pg = pgls(y, X, C, ["intercept", pname])
                entry["pgls_" + pname] = {
                    "beta": pg["coef"][pname]["beta"],
                    "t": pg["coef"][pname]["t"],
                    "p_two_sided": pg["coef"][pname]["p_two_sided"],
                    "r2_gls": pg["r2_gls"], "aicc": pg.get("aicc")}
            sub[resp] = entry
        phy[rootmode] = sub
    # rooting-free controls
    D = patristic(tree, present)
    rf = {}
    for resp in RESPONSES:
        y = np.array([d[resp] for d in rows], float)
        rf[resp] = {"mantel_patristic_vs_traitdiff": mantel(D, y),
                    "morans_I_inv_patristic": morans_i(D, y)}
    phy["rooting_free"] = rf
    phy["patristic_matrix"] = {"tips": present, "D": D.tolist()}
    rec["phylogenetic_control"] = phy

    # ---- (6) leave-one-species-out on the headline test
    loo = []
    for drop in range(k):
        keep = [i for i in range(k) if i != drop]
        y = np.array([rows[i]["piN_piS"] for i in keep], float)
        g = np.array([rows[i]["category"] == "pathogenic" for i in keep], bool)
        if 0 < g.sum() < len(keep):
            m = mwu_exact(y, g)
            loo.append({"dropped": tips[drop], "k": len(keep),
                        "n1": m["n1"], "n2": m["n2"], "U1": m["U1"],
                        "p_two_sided": m["p_two_sided"],
                        "p_floor": m["p_floor_two_sided"],
                        "median_path": m["median_group1"],
                        "median_env": m["median_group2"],
                        "direction_median": ("pathogenic_lower"
                                             if m["median_group1"] < m["median_group2"]
                                             else "pathogenic_higher"),
                        "mean_path": m["mean_group1"],
                        "mean_env": m["mean_group2"],
                        "direction_mean": ("pathogenic_lower"
                                           if m["mean_group1"] < m["mean_group2"]
                                           else "pathogenic_higher")})
    rec["leave_one_species_out_piN_piS_category"] = loo

    # ---- (7) stratified deconfounding: clade within disease, disease within host
    strat = []
    strata = [
        ("within_pathogenic__clade", lambda d: d["category"] == "pathogenic",
         lambda d: d["nw_hps_clade"],
         "among PATHOGENIC species only: New World HPS clade vs Old World. "
         "Isolates clade from disease."),
        ("within_rodent__disease", lambda d: d["is_rodent"],
         lambda d: d["category"] == "pathogenic",
         "among RODENT-borne species only: pathogenic vs environmental. "
         "Isolates disease from host class."),
        ("within_nonNWclade__disease", lambda d: not d["nw_hps_clade"],
         lambda d: d["category"] == "pathogenic",
         "among species OUTSIDE the New World HPS clade: pathogenic vs "
         "environmental.  Isolates disease from clade."),
        ("within_rodent_nonNWclade__disease",
         lambda d: d["is_rodent"] and not d["nw_hps_clade"],
         lambda d: d["category"] == "pathogenic",
         "rodent-borne AND outside the New World HPS clade: the only stratum in "
         "which disease is free of both confounds."),
    ]
    for sname, keepf, groupf, desc in strata:
        sub = [d for d in rows if keepf(d)]
        entry = {"stratum": sname, "description": desc, "k": len(sub),
                 "members": [d["abbrev"] for d in sub]}
        if len(sub) >= 3:
            g = np.array([groupf(d) for d in sub], bool)
            entry["n_group1"] = int(g.sum())
            entry["n_group2"] = len(sub) - int(g.sum())
            if 0 < g.sum() < len(sub):
                for resp in RESPONSES:
                    y = np.array([d[resp] for d in sub], float)
                    entry[resp] = mwu_exact(y, g)
            else:
                entry["testable"] = False
                entry["reason"] = "one group is empty in this stratum"
        else:
            entry["testable"] = False
            entry["reason"] = "fewer than 3 species in the stratum"
        strat.append(entry)
    rec["stratified_deconfounding"] = strat

    out_json[dsname] = rec
    return rec


def mde_block(rows, outdir):
    """Minimum detectable effect, three independent ways."""
    k = len(rows)
    y = np.array([d["piN_piS"] for d in rows], float)
    g = np.array([d["category"] == "pathogenic" for d in rows], bool)
    tips = [d["abbrev"] for d in rows]
    k1, k2 = int(g.sum()), k - int(g.sum())
    yp, ye = y[g], y[~g]
    out = {"k_pathogenic": k1, "k_environmental": k2,
           "observed_values": {t: float(v) for t, v in zip(tips, y)}}

    # (i) exact-test p floor
    n_assign = math.comb(k, k2)
    out["exact_test_p_floor"] = {
        "n_label_assignments": n_assign,
        "min_attainable_two_sided_p": 2.0 / n_assign,
        "requires": "perfect, non-overlapping separation of the two groups",
        "alpha": TWO_SIDED_ALPHA,
        "any_p_below_alpha_attainable": bool(2.0 / n_assign < TWO_SIDED_ALPHA)}

    # (ii) exact separation requirement, holding the pathogenic values fixed
    need_above = float(yp.max())
    need_below = float(yp.min())
    out["separation_requirement"] = {
        "direction_pathogenic_lower": {
            "condition": "min(environmental) > max(pathogenic) = %.6g" % need_above,
            "required_fold_change_per_species": {
                t: (None if v <= 0 else float(need_above / v))
                for i, (t, v) in enumerate(zip(tips, y)) if not g[i]},
        },
        "direction_pathogenic_higher": {
            "condition": "max(environmental) < min(pathogenic) = %.6g" % need_below,
            "attainable": bool(need_below > 0),
            "note": ("impossible: the global minimum piN/piS is a pathogenic "
                     "species (%s = %.6g)" % (tips[int(np.argmin(y))], y.min()))
            if need_below <= min(ye) else ""}}

    # (iii) parametric two-sample t MDE on the raw scale
    s2p = (((k1 - 1) * yp.var(ddof=1) + (k2 - 1) * ye.var(ddof=1))
           / (k1 + k2 - 2))
    sp = math.sqrt(s2p)
    dof = k1 + k2 - 2
    # t critical value by bisection on the Student CDF
    lo, hi = 0.0, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if student_p(mid, dof) > TWO_SIDED_ALPHA:
            lo = mid
        else:
            hi = mid
    tcrit = 0.5 * (lo + hi)
    mde = tcrit * sp * math.sqrt(1.0 / k1 + 1.0 / k2)
    obs_diff = float(ye.mean() - yp.mean())
    out["parametric_t_mde"] = {
        "pooled_sd": sp, "dof": dof, "t_crit_two_sided_0.05": tcrit,
        "mde_difference_in_means_raw": mde,
        "mde_as_fold_of_pathogenic_mean": mde / float(yp.mean()),
        "observed_difference_env_minus_path": obs_diff,
        "observed_over_mde": obs_diff / mde,
        "pathogenic_mean": float(yp.mean()), "environmental_mean": float(ye.mean()),
        "pathogenic_sd": float(yp.std(ddof=1)),
        "environmental_sd": float(ye.std(ddof=1)),
        "cohens_d": obs_diff / sp}

    # (iv) simulation power curve: lognormal model fitted to the nonzero
    #      pathogenic values, environmental group multiplied by a fold change f
    nz = yp[yp > 0]
    mu = float(np.log10(nz).mean())
    sd = float(np.log10(nz).std(ddof=1))
    rng = np.random.default_rng([SEED, 97])
    curve = []
    folds = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0]
    nsim = 5000
    for f in folds:
        hit = 0
        for _ in range(nsim):
            a = 10 ** (rng.normal(mu, sd, k1))
            b = 10 ** (rng.normal(mu + math.log10(f), sd, k2))
            vals = np.concatenate([a, b])
            gg = np.array([True] * k1 + [False] * k2)
            if mwu_exact(vals, gg)["p_two_sided"] <= TWO_SIDED_ALPHA:
                hit += 1
        curve.append({"fold_change": f, "power": hit / nsim, "n_sim": nsim})
    out["simulation_power_curve"] = {
        "model": ("log10(piN/piS) ~ Normal(mu, sd) fitted to the %d NONZERO "
                  "pathogenic values; environmental group shifted by log10(f)" % len(nz)),
        "mu_log10": mu, "sd_log10": sd, "k1": k1, "k2": k2,
        "n_sim_per_point": nsim, "seed": SEED, "alpha": TWO_SIDED_ALPHA,
        "curve": curve,
        "fold_change_for_80pct_power": next(
            (c["fold_change"] for c in curve if c["power"] >= 0.80), None)}
    return out


def main():
    root, labels_tsv, tree_nwk, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    labels = load_labels(labels_tsv)
    tree = read_tree(tree_nwk)
    with open(os.path.join(outdir, "effort_covariates.json")) as fh:
        effort = json.load(fh)
    table = load_table(os.path.join(root, "results", "pinps_per_species.tsv"))
    design, nw_clade = build_design(root, labels, tree, effort, table)

    datasets = [
        ("A_primary_clust95d_n10", "clust95d", "10", True,
         "PRIMARY: single-population sets, all species subsampled to n=10 "
         "(equal study effort by construction), 1000 draws, seed 20261012"),
        ("B_clust95d_full", "clust95d", "full", True,
         "same 10 species at their full single-population n (11-291): the only "
         "level at which the effort regression has variation in n actually used"),
        ("C_clust95d_n12", "clust95d", "12", True,
         "equal n=12 (loses CHOV); k=9"),
        ("D_clust95_full", "clust95", "full", True,
         "single-linkage clusters instead of diameter-pruned; loses PUUV and TULV "
         "to the lineage-pool flag, so NO rodent-borne environmental remains"),
        ("E_ss_all_n10_SATURATED", "ss_all", "full", False,
         "SENSITIVITY ONLY - whole ICTV species, n>=10, piS SATURATED "
         "(piS_jc ~1); included because it is the only level with k=12 and "
         "4 non-rodent environmentals"),
        ("F_ss_uniq_n10_SATURATED", "ss_uniq", "full", False,
         "SENSITIVITY ONLY - as E but one sequence per haplotype; also saturated"),
    ]

    out_tests, out_json = [], {}
    first_rows = None
    for name, ss, lv, usable, note in datasets:
        rows = dataset(table, design, ss, lv, usable)
        if len(rows) < 4:
            out_json[name] = {"dataset": name, "note": note, "k": len(rows),
                              "skipped": "fewer than 4 species"}
            continue
        print("== %s: k=%d (%s)" % (name, len(rows),
                                    ",".join(d["abbrev"] for d in rows)), flush=True)
        run_battery(rows, name, note, tree, out_tests, out_json)
        if first_rows is None:
            first_rows = rows

    out_json["_meta"] = {
        "seed": SEED,
        "nw_hps_clade_from_tree": nw_clade,
        "nw_hps_clade_definition":
            "all tips of MRCA(SNV, BAYV) in S22_ml_clean.nwk",
        "tree": tree_nwk,
        "labels": labels_tsv,
        "effort_source": "results/stats_confound-direct/effort_covariates.json "
                         "(NCBI esearch, retrieved %s)" % effort["retrieved_utc"],
        "responses": RESPONSES,
        "binary_predictors": [b[0] for b in BIN_PREDICTORS],
        "continuous_predictors": [c[0] for c in CONT_PREDICTORS],
    }
    out_json["_mde"] = mde_block(first_rows, outdir)
    out_json["_design_matrix"] = first_rows

    with open(os.path.join(outdir, "stats_confound_direct.json"), "w") as fh:
        json.dump(out_json, fh, indent=1, sort_keys=True, default=str)

    # tidy TSVs
    keys = []
    for r in out_tests:
        for kk in r:
            if kk not in keys:
                keys.append(kk)
    with open(os.path.join(outdir, "tests_all.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in out_tests:
            w.writerow(r)

    ve_rows = []
    for name, rec in out_json.items():
        if isinstance(rec, dict) and "variance_explained" in rec:
            ve_rows.extend(rec["variance_explained"])
    with open(os.path.join(outdir, "variance_explained.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ve_rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in ve_rows:
            w.writerow(r)

    rs_rows = []
    for name, rec in out_json.items():
        if isinstance(rec, dict) and "residual_tests" in rec:
            rs_rows.extend(rec["residual_tests"])
    with open(os.path.join(outdir, "residual_tests.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rs_rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rs_rows:
            w.writerow(r)

    dm_keys = list(first_rows[0].keys())
    with open(os.path.join(outdir, "design_matrix.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=dm_keys, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in first_rows:
            w.writerow(r)
    print("wrote", outdir, flush=True)


if __name__ == "__main__":
    main()
