"""Phylogenetic statistical treatment of within-species piN/piS (nucleoprotein).

Fits phylogenetic generalised least squares (PGLS) and Garland-style phylogenetic
ANOVA of per-species piN/piS against (a) disease category, (b) reservoir host
class, (c) both, and (d) the New World HPS clade, on the 22-taxon ML tree pruned
to the species that have a usable within-species value.

Everything is hand-rolled on numpy/scipy (no ape, no phytools, no ete3).

Usage:
  python3 -I 21_phylo_stats.py <ANALYSIS_ROOT> <LABELS_TSV> <TREE_NWK> <OUTDIR>

Deterministic: master seed 20261012.
"""
import csv
import itertools
import json
import math
import os
import sys
import zlib

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "phylo_tree", os.path.join(os.path.dirname(os.path.abspath(__file__)), "20_phylo_tree.py")
)
pt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pt)

SEED = 20261012
NSIM_ANOVA = 10000
NSIM_POWER = 4000

NW_HPS_CLADE = {"SNV", "ANDV", "LANV", "CHOV", "BCCV", "BAYV"}


# ============================================================ PGLS machinery
def vmat(C, lam):
    """Pagel's lambda transform of a BM VCV: off-diagonals scaled by lambda."""
    V = lam * C.copy()
    np.fill_diagonal(V, np.diag(C))
    return V


def gls(y, X, V):
    """ML GLS fit. Returns dict with beta, se, t, p, logL, sigma2, resid."""
    n = len(y)
    p = X.shape[1]
    L = np.linalg.cholesky(V)
    Xw = np.linalg.solve(L, X)
    yw = np.linalg.solve(L, y)
    XtX = Xw.T @ Xw
    XtXi = np.linalg.inv(XtX)
    beta = XtXi @ (Xw.T @ yw)
    r = yw - Xw @ beta
    rss = float(r @ r)
    sigma2_ml = rss / n
    sigma2_reml = rss / (n - p) if n > p else np.nan
    logdetV = 2.0 * float(np.sum(np.log(np.diag(L))))
    logL = -0.5 * (n * math.log(2 * math.pi) + n * math.log(sigma2_ml) + logdetV + n)
    cov = sigma2_reml * XtXi
    se = np.sqrt(np.diag(cov))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = beta / se
    df = n - p
    pv = 2 * stats.t.sf(np.abs(t), df) if df > 0 else np.full(p, np.nan)
    return dict(
        beta=beta, se=se, t=t, p=pv, logL=logL, sigma2_ml=sigma2_ml,
        sigma2_reml=sigma2_reml, rss=rss, df=df, npar_mean=p,
    )


def fit_lambda(y, X, C, bounds=(0.0, 1.0)):
    """ML estimate of Pagel's lambda jointly with the GLS fit."""
    def neg(lam):
        return -gls(y, X, vmat(C, lam))["logL"]

    grid = np.linspace(bounds[0], bounds[1], 101)
    vals = [neg(g) for g in grid]
    g0 = grid[int(np.argmin(vals))]
    lo = max(bounds[0], g0 - 0.02)
    hi = min(bounds[1], g0 + 0.02)
    res = minimize_scalar(neg, bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-6})
    lam = float(res.x) if res.fun <= min(vals) else float(g0)
    fit = gls(y, X, vmat(C, lam))
    fit["lambda"] = lam
    return fit


def aic_set(logL, k, n):
    aic = -2 * logL + 2 * k
    if n - k - 1 > 0:
        aicc = aic + 2 * k * (k + 1) / (n - k - 1)
    else:
        aicc = float("nan")
    return aic, aicc


def blomberg_K(y, C):
    """Blomberg's K: observed MSE0/MSE ratio relative to the BM expectation."""
    n = len(y)
    one = np.ones((n, 1))
    Vi = np.linalg.inv(C)
    a = float((one.T @ Vi @ y) / (one.T @ Vi @ one))
    r = y - a
    mse = float(r @ Vi @ r) / (n - 1)
    mse0 = float(r @ r) / (n - 1)
    expected = (np.trace(C) - n / float(np.sum(Vi))) / (n - 1)
    return (mse0 / mse) / expected


# ============================================================ data loading
def dseed(*parts):
    """Deterministic integer seed component from strings (hash() is randomised)."""
    return [int(zlib.crc32(str(x).encode())) for x in parts]


def load_labels(path):
    lab = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            lab[row["abbrev"]] = row
    return lab


def load_table(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pick(rows, seqset, level, require_usable):
    out = []
    for r in rows:
        if r["seqset"] != seqset or r["n_level"] != level:
            continue
        if require_usable and r["usable_for_group_test"] != "True":
            continue
        out.append(r)
    return out


# ============================================================ design matrices
def design(species, labelling, labels):
    """Return (X, colnames, group_vector) for a labelling."""
    n = len(species)
    if labelling == "intercept":
        return np.ones((n, 1)), ["intercept"], None
    if labelling == "disease":
        g = np.array([1.0 if labels[s]["category"] == "environmental" else 0.0
                      for s in species])
        return np.column_stack([np.ones(n), g]), ["intercept", "environmental"], g
    if labelling == "disease_tulv_moved":
        g = np.array([1.0 if (labels[s]["category"] == "environmental" and s != "TULV")
                      else 0.0 for s in species])
        return np.column_stack([np.ones(n), g]), ["intercept", "environmental_tulvmoved"], g
    if labelling == "host2":
        g = np.array([0.0 if labels[s]["is_rodent"] == "True" else 1.0 for s in species])
        return np.column_stack([np.ones(n), g]), ["intercept", "non_rodent"], g
    if labelling == "host3":
        hc = [labels[s]["host_class"] for s in species]
        lev = sorted(set(hc))
        base = "rodent" if "rodent" in lev else lev[0]
        others = [l for l in lev if l != base]
        cols = [np.ones(n)]
        names = ["intercept"]
        for o in others:
            cols.append(np.array([1.0 if h == o else 0.0 for h in hc]))
            names.append("host_" + o)
        g = np.array([lev.index(h) for h in hc], dtype=float)
        return np.column_stack(cols), names, g
    if labelling == "cladeNW":
        g = np.array([1.0 if s in NW_HPS_CLADE else 0.0 for s in species])
        return np.column_stack([np.ones(n), g]), ["intercept", "NW_HPS_clade"], g
    if labelling == "disease+host2":
        d = np.array([1.0 if labels[s]["category"] == "environmental" else 0.0
                      for s in species])
        h = np.array([0.0 if labels[s]["is_rodent"] == "True" else 1.0 for s in species])
        return (np.column_stack([np.ones(n), d, h]),
                ["intercept", "environmental", "non_rodent"], None)
    raise ValueError(labelling)


# ============================================================ tests
def anova_F(y, g):
    """One-way ANOVA F for grouping vector g (floats as group codes)."""
    groups = [y[g == v] for v in np.unique(g)]
    groups = [a for a in groups if len(a) > 0]
    k = len(groups)
    n = len(y)
    if k < 2 or n - k <= 0:
        return float("nan")
    gm = y.mean()
    ssb = sum(len(a) * (a.mean() - gm) ** 2 for a in groups)
    ssw = sum(float(((a - a.mean()) ** 2).sum()) for a in groups)
    if ssw == 0:
        return float("inf")
    return (ssb / (k - 1)) / (ssw / (n - k))


def phylo_anova(y, g, C, rng, nsim=NSIM_ANOVA):
    """Garland et al. (1993) phylogenetic ANOVA: conventional F, null from BM sims."""
    n = len(y)
    one = np.ones((n, 1))
    f0 = gls(y, one, C)
    mu = float(f0["beta"][0])
    sigma2 = f0["sigma2_reml"]
    F_obs = anova_F(y, g)
    k = len(np.unique(g))
    p_param = 1.0 - stats.f.cdf(F_obs, k - 1, n - k) if np.isfinite(F_obs) else np.nan
    L = np.linalg.cholesky(sigma2 * C)
    Z = rng.standard_normal((n, nsim))
    Y = mu + L @ Z
    Fs = np.array([anova_F(Y[:, i], g) for i in range(nsim)])
    p_phylo = float((np.sum(Fs >= F_obs) + 1) / (nsim + 1))
    return dict(F_obs=float(F_obs), p_ols=float(p_param), p_phylo=p_phylo,
                nsim=nsim, mu=mu, sigma2=float(sigma2))


def mannwhitney_exact(a, b):
    """Exact two-sided Mann-Whitney U (no ties assumed; ties handled by midranks
    for the statistic but the exact p uses the no-tie permutation null)."""
    try:
        res = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
    except Exception:
        res = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(res.statistic), float(res.pvalue)


def exact_label_randomisation(y, C, species, g_obs, labels, tree_root,
                              restrict_parsimony=False):
    """Enumerate every assignment of the observed group sizes to the tips, refit
    PGLS (BM) for each, and compare |t| of the group coefficient.

    With restrict_parsimony=True, only assignments whose parsimony score on the
    tree equals the observed labelling's score are retained — a
    Maddison-concentrated-changes-style randomisation that preserves how
    phylogenetically clustered the labels are.
    """
    n = len(y)
    m = int(g_obs.sum())
    idx_obs = tuple(sorted(np.where(g_obs == 1)[0].tolist()))
    states_obs = {species[i]: int(g_obs[i]) for i in range(n)}
    pars_obs = pt.parsimony_score(tree_root, states_obs)

    L = np.linalg.cholesky(C)
    yw = np.linalg.solve(L, y)
    onew = np.linalg.solve(L, np.ones((n, 1)))

    def tstat(idx):
        g = np.zeros(n)
        g[list(idx)] = 1.0
        gw = np.linalg.solve(L, g.reshape(-1, 1))
        Xw = np.hstack([onew, gw])
        XtXi = np.linalg.inv(Xw.T @ Xw)
        beta = XtXi @ (Xw.T @ yw)
        r = yw - Xw @ beta
        s2 = float(r @ r) / (n - 2)
        se = math.sqrt(s2 * XtXi[1, 1])
        return abs(float(beta[1]) / se) if se > 0 else float("inf"), float(beta[1])

    t_obs, b_obs = tstat(idx_obs)
    kept, ts, diffs = [], [], []
    for idx in itertools.combinations(range(n), m):
        st = {species[i]: (1 if i in idx else 0) for i in range(n)}
        ps = pt.parsimony_score(tree_root, st)
        if restrict_parsimony and ps != pars_obs:
            continue
        ti, bi = tstat(idx)
        kept.append(idx)
        ts.append(ti)
        sel = np.zeros(n, bool)
        sel[list(idx)] = True
        diffs.append(float(y[sel].mean() - y[~sel].mean()))
    ts = np.array(ts)
    diffs = np.array(diffs)
    d_obs = float(y[g_obs == 1].mean() - y[g_obs == 0].mean())
    return dict(
        n_assignments=len(kept),
        parsimony_obs=pars_obs,
        restricted=bool(restrict_parsimony),
        t_obs=float(t_obs), beta_obs=float(b_obs),
        p_t=float(np.sum(ts >= t_obs - 1e-12) / len(ts)),
        diff_obs=d_obs,
        p_diff=float(np.sum(np.abs(diffs) >= abs(d_obs) - 1e-12) / len(diffs)),
        min_attainable_p=float(1.0 / len(ts)),
    )


def power_curve(C, g, rng, sigma2, d_grid, nsim=NSIM_POWER, alpha=0.05):
    """Power of (a) the exact label-randomisation test on the BM-PGLS t statistic
    and (b) the exact two-sided Mann-Whitney test, as a function of the true
    group difference d, with residuals simulated as BM on the tree."""
    n = len(g)
    m = int(g.sum())
    L = np.linalg.cholesky(sigma2 * C)
    Lc = np.linalg.cholesky(C)
    onew = np.linalg.solve(Lc, np.ones((n, 1)))

    combos = list(itertools.combinations(range(n), m))
    projs = []
    for idx in combos:
        gg = np.zeros(n)
        gg[list(idx)] = 1.0
        gw = np.linalg.solve(Lc, gg.reshape(-1, 1))
        Xw = np.hstack([onew, gw])
        XtXi = np.linalg.inv(Xw.T @ Xw)
        P = XtXi @ Xw.T
        H = Xw @ P
        projs.append((P, np.eye(n) - H, XtXi[1, 1]))
    obs_pos = combos.index(tuple(sorted(np.where(g == 1)[0].tolist())))

    env = g == 1
    out = []
    for d in d_grid:
        Y = (d * g).reshape(-1, 1) + L @ rng.standard_normal((n, nsim))
        Yw = np.linalg.solve(Lc, Y)
        T = np.empty((len(combos), nsim))
        for j, (P, M, v) in enumerate(projs):
            b1 = (P @ Yw)[1, :]
            R = M @ Yw
            s2 = np.einsum("ij,ij->j", R, R) / (n - 2)
            T[j, :] = np.abs(b1) / np.sqrt(s2 * v)
        t_obs = T[obs_pos, :]
        p_rand = (T >= t_obs - 1e-12).sum(axis=0) / len(combos)
        rej_rand = float(np.mean(p_rand <= alpha + 1e-12))
        # exact Mann-Whitney: with these group sizes only perfect separation can
        # reach alpha=0.05, so count perfect separations
        a = Y[env, :]
        b = Y[~env, :]
        perfect = np.mean((a.max(axis=0) < b.min(axis=0)) | (a.min(axis=0) > b.max(axis=0)))
        out.append(dict(d=float(d), power_randomisation=rej_rand,
                        power_mannwhitney=float(perfect)))
    return out


# ============================================================ main
def run_dataset(tag, rows, labels, tree_full, outdir, notes):
    species = sorted(r["abbrev"] for r in rows)
    byab = {r["abbrev"]: r for r in rows}
    pruned = pt.prune(tree_full, set(species))
    pruned.parent = None
    adj, names = pt.adjacency(pruned)

    # ---- rootings
    u, v, Lb, frac, span = pt.midpoint_edge(adj, names, species)
    rootings = {}
    rootings["midpoint"] = pt.vcv_rooted_at_edge(adj, names, species, u, v, Lb, frac)
    # root at the top node of the newick as given (IQ-TREE's arbitrary placement)
    name_to_id = {nm: k for k, nm in names.items() if nm}
    rootings["as_given"] = pt.vcv_rooted_at_node(adj, names, species, 0)
    # outgroup rooting on the non-rodent-borne clade, if it is monophyletic here
    nonrod = [s for s in species if labels[s]["is_rodent"] != "True"]
    og = pt.find_edge_subtending(pruned, nonrod) if nonrod else None
    if og is not None and og.parent is not None:
        # rebuild ids by walking pruned tree in the same order as adjacency()
        nodes = []
        stack = [pruned]
        while stack:
            nn = stack.pop()
            nodes.append(nn)
            stack.extend(nn.children)
        ids = {id(nn): i for i, nn in enumerate(nodes)}
        nid = ids[id(og)]
        pid = ids[id(og.parent)]
        rootings["outgroup_nonrodent"] = pt.vcv_rooted_at_edge(
            adj, names, species, pid, nid, og.blen, 0.5
        )

    # ---- responses
    resp = {}
    pp = np.array([float(byab[s]["piN_piS"]) for s in species])
    resp["piN_piS"] = pp
    nz = pp[pp > 0]
    off = float(nz.min()) / 2.0
    resp["log10_piN_piS_off"] = np.log10(pp + off)
    resp["piN"] = np.array([float(byab[s]["piN"]) for s in species])
    resp["piS"] = np.array([float(byab[s]["piS"]) for s in species])
    resp["rank_piN_piS"] = stats.rankdata(pp)

    n_used = {s: int(byab[s]["n_sequences"]) for s in species}

    result = dict(
        tag=tag, species=species, k_species=len(species),
        n_sequences_used=n_used,
        category={s: labels[s]["category"] for s in species},
        host_class={s: labels[s]["host_class"] for s in species},
        is_rodent={s: labels[s]["is_rodent"] for s in species},
        NW_HPS_clade={s: (s in NW_HPS_CLADE) for s in species},
        piN_piS={s: float(byab[s]["piN_piS"]) for s in species},
        piN={s: float(byab[s]["piN"]) for s in species},
        piS={s: float(byab[s]["piS"]) for s in species},
        log_offset=off,
        midpoint_span=float(span),
        rootings=sorted(rootings),
        notes=notes,
        fits=[], signal=[], anova=[], randomisation=[], power=None,
        crosstab={}, saturation_flags={s: byab[s]["flag_piS_saturated"] for s in species},
    )

    # collinearity between the disease labelling and the host labelling
    _, _, gd = design(species, "disease", labels)
    _, _, gh = design(species, "host2", labels)
    _, _, gc = design(species, "cladeNW", labels)
    result["collinearity"] = dict(
        pearson_disease_host2=float(np.corrcoef(gd, gh)[0, 1]),
        pearson_disease_cladeNW=float(np.corrcoef(gd, gc)[0, 1]),
        n_nonrodent=int(gh.sum()), n_environmental=int(gd.sum()),
        n_environmental_rodent=int(np.sum((gd == 1) & (gh == 0))),
        n_pathogenic_nonrodent=int(np.sum((gd == 0) & (gh == 1))),
        discriminating_species=[species[i] for i in range(len(species))
                                if gd[i] != gh[i]],
        parsimony_disease=pt.parsimony_score(
            pruned, {species[i]: int(gd[i]) for i in range(len(species))}),
        parsimony_host2=pt.parsimony_score(
            pruned, {species[i]: int(gh[i]) for i in range(len(species))}),
        parsimony_cladeNW=pt.parsimony_score(
            pruned, {species[i]: int(gc[i]) for i in range(len(species))}),
    )
    result["pruned_newick"] = pt.to_newick(pruned)

    # confound crosstab
    ct = {}
    for s in species:
        key = labels[s]["category"] + " x " + labels[s]["host_class"]
        ct[key] = ct.get(key, 0) + 1
    result["crosstab"] = ct

    labellings = ["intercept", "disease", "disease_tulv_moved", "host2", "host3",
                  "disease+host2", "cladeNW"]

    for rootname, C in sorted(rootings.items()):
        for yname, y in resp.items():
            # phylogenetic signal of the response
            one = np.ones((len(y), 1))
            fB = gls(y, one, C)
            fO = gls(y, one, vmat(C, 0.0))
            fL = fit_lambda(y, one, C)
            lrt = 2 * (fL["logL"] - fO["logL"])
            result["signal"].append(dict(
                rooting=rootname, response=yname,
                lambda_ML=fL["lambda"],
                logL_BM=fB["logL"], logL_lambda=fL["logL"], logL_star=fO["logL"],
                LRT_lambda_vs_0=float(lrt),
                p_LRT_lambda_vs_0=float(stats.chi2.sf(max(lrt, 0.0), 1)),
                blomberg_K=float(blomberg_K(y, C)),
                aicc_BM=aic_set(fB["logL"], 2, len(y))[1],
                aicc_star=aic_set(fO["logL"], 2, len(y))[1],
                aicc_lambda=aic_set(fL["logL"], 3, len(y))[1],
            ))
            for lab in labellings:
                X, cn, g = design(species, lab, labels)
                if np.linalg.matrix_rank(X) < X.shape[1]:
                    result["fits"].append(dict(
                        rooting=rootname, response=yname, model=lab,
                        status="RANK DEFICIENT - not identifiable",
                        colnames=cn))
                    continue
                for model, V, extra_par in (
                    ("BM", C, 0), ("lambda_ML", None, 1), ("star_OLS", vmat(C, 0.0), 0)
                ):
                    if model == "lambda_ML":
                        f = fit_lambda(y, X, C)
                        lam = f["lambda"]
                    else:
                        f = gls(y, X, V)
                        lam = 1.0 if model == "BM" else 0.0
                    k = X.shape[1] + 1 + extra_par
                    aic, aicc = aic_set(f["logL"], k, len(y))
                    rec = dict(rooting=rootname, response=yname, model=lab,
                               error_model=model, lam=float(lam),
                               colnames=cn, k_params=k,
                               logL=f["logL"], AIC=aic, AICc=aicc,
                               beta=[float(x) for x in f["beta"]],
                               se=[float(x) for x in f["se"]],
                               t=[float(x) for x in f["t"]],
                               p=[float(x) for x in f["p"]],
                               df=f["df"], sigma2_reml=float(f["sigma2_reml"]))
                    # LRT against intercept-only with same error model
                    if lab != "intercept":
                        if model == "lambda_ML":
                            f0 = fit_lambda(y, one, C)
                        else:
                            f0 = gls(y, one, V)
                        d = 2 * (f["logL"] - f0["logL"])
                        ddf = X.shape[1] - 1
                        rec["LRT_vs_intercept"] = float(d)
                        rec["LRT_df"] = ddf
                        rec["p_LRT"] = float(stats.chi2.sf(max(d, 0.0), ddf))
                    result["fits"].append(rec)

    # ---- Garland phylogenetic ANOVA, under every rooting
    for rootname, Cp in sorted(rootings.items()):
        for yname, y in resp.items():
            for lab in ["disease", "disease_tulv_moved", "host2", "host3", "cladeNW"]:
                X, cn, g = design(species, lab, labels)
                if g is None or len(np.unique(g)) < 2:
                    continue
                rng = np.random.default_rng(
                    [SEED, 1] + dseed(tag, rootname, yname, lab))
                r = phylo_anova(y, g, Cp, rng)
                r.update(rooting=rootname, response=yname, model=lab,
                         group_sizes=[int(np.sum(g == v)) for v in np.unique(g)])
                result["anova"].append(r)

    # ---- exact label randomisation on the PGLS t statistic, under every rooting
    for rootname, Cp in sorted(rootings.items()):
        for yname, y in resp.items():
            for lab in ["disease", "disease_tulv_moved", "cladeNW", "host2"]:
                X, cn, g = design(species, lab, labels)
                if g is None or g.sum() in (0, len(g)):
                    continue
                for restrict in (False, True):
                    r = exact_label_randomisation(y, Cp, species, g, labels, pruned,
                                                  restrict_parsimony=restrict)
                    r.update(rooting=rootname, response=yname, model=lab)
                    # non-phylogenetic exact Mann-Whitney for reference
                    a = y[g == 1]
                    b = y[g == 0]
                    U, pU = mannwhitney_exact(a, b)
                    r["mannwhitney_U"] = U
                    r["mannwhitney_p_exact"] = pU
                    r["group_sizes"] = [int(len(a)), int(len(b))]
                    result["randomisation"].append(r)

    # ---- power / minimum detectable effect, primary response, every rooting
    result["power"] = []
    y = resp["piN_piS"]
    X, cn, g = design(species, "disease", labels)
    obs_diff = float(y[g == 1].mean() - y[g == 0].mean())
    for rootname, Cp in sorted(rootings.items()):
        f0 = gls(y, np.ones((len(y), 1)), Cp)
        sigma2 = f0["sigma2_reml"]
        rng = np.random.default_rng([SEED, 2] + dseed(tag, rootname))
        d_grid = [0.0, 0.005, 0.01, obs_diff, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10,
                  0.15, 0.20, 0.30, 0.50]
        d_grid = sorted(set(round(d, 6) for d in d_grid))
        result["power"].append(dict(
            rooting=rootname, response="piN_piS",
            sigma2_reml_BM=float(sigma2),
            implied_tip_sd=float(np.sqrt(sigma2 * np.mean(np.diag(Cp)))),
            observed_diff_env_minus_path=obs_diff,
            nsim=NSIM_POWER, alpha=0.05, seed=SEED,
            curve=power_curve(Cp, g, rng, sigma2, d_grid),
        ))
    return result


def main():
    root, labels_tsv, tree_nwk, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    labels = load_labels(labels_tsv)
    rows = load_table(os.path.join(root, "results", "pinps_per_species.tsv"))
    tree_full = pt.parse_newick(open(tree_nwk).read())

    datasets = []

    # PRIMARY: clust95d at equal depth n = 10 (all ten usable species)
    r1 = pick(rows, "clust95d", "10", True)
    datasets.append(("clust95d_n10", r1,
                     "PRIMARY. Single-population sets (<=5% diameter), every species "
                     "subsampled to n = 10 sequences, 1000 draws, seed 20261012. "
                     "Study effort equalised."))

    # SENSITIVITY A: clust95d full depth
    r2 = pick(rows, "clust95d", "full", True)
    datasets.append(("clust95d_full", r2,
                     "SENSITIVITY. Same species, full sequence depth (n = 12-291), "
                     "so study effort is NOT equalised."))

    # SENSITIVITY B: ss_all, every species with n >= 10 (premise fails: piS saturated)
    r3 = [r for r in pick(rows, "ss_all", "full", False)
          if int(r["n_sequences"]) >= 10]
    datasets.append(("ss_all_n_ge_10", r3,
                     "SENSITIVITY, PREMISE INVALID. Whole sensu-stricto species sets "
                     "with n >= 10 (12 species). 11 of the 12 have "
                     "flag_piS_saturated = True (CHOV is the exception), i.e. piS is "
                     "saturated and piN/piS is not interpretable. Included only "
                     "because it is the one set with more than one species in a "
                     "non-rodent host-class level."))

    # SENSITIVITY C: clust95 (single-linkage, not diameter-pruned)
    r4 = pick(rows, "clust95", "full", True)
    datasets.append(("clust95_full", r4,
                     "SENSITIVITY. Single-linkage 95% clusters without the diameter "
                     "prune; PUUV and TULV drop out on flag_lineage_pool."))

    out = dict(seed=SEED, tree=tree_nwk, labels=labels_tsv,
               table=os.path.join(root, "results", "pinps_per_species.tsv"),
               datasets=[])
    for tag, rws, note in datasets:
        if len(rws) < 4:
            out["datasets"].append(dict(tag=tag, skipped=True, k_species=len(rws),
                                        notes=note))
            continue
        out["datasets"].append(run_dataset(tag, rws, labels, tree_full, outdir, note))

    with open(os.path.join(outdir, "phylo_stats.json"), "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    # tidy TSV of every model fit
    with open(os.path.join(outdir, "pgls_model_fits.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dataset", "k_species", "rooting", "response", "model",
                    "error_model", "lambda", "k_params", "logL", "AIC", "AICc",
                    "beta_group", "se_group", "t_group", "p_group", "df",
                    "LRT_vs_intercept", "LRT_df", "p_LRT", "status"])
        for ds in out["datasets"]:
            if ds.get("skipped"):
                continue
            for f in ds["fits"]:
                if f.get("status"):
                    w.writerow([ds["tag"], ds["k_species"], f["rooting"], f["response"],
                                f["model"], "", "", "", "", "", "", "", "", "", "",
                                "", "", "", "", f["status"]])
                    continue
                gi = 1 if len(f["beta"]) > 1 else 0
                w.writerow([ds["tag"], ds["k_species"], f["rooting"], f["response"],
                            f["model"], f["error_model"], "%.4f" % f["lam"],
                            f["k_params"], "%.4f" % f["logL"], "%.3f" % f["AIC"],
                            "%.3f" % f["AICc"],
                            "%.6g" % f["beta"][gi], "%.6g" % f["se"][gi],
                            "%.4f" % f["t"][gi], "%.4g" % f["p"][gi], f["df"],
                            "%.4f" % f.get("LRT_vs_intercept", float("nan")),
                            f.get("LRT_df", ""),
                            "%.4g" % f.get("p_LRT", float("nan")), "OK"])
    # tidy TSV of the per-species input data actually used
    with open(os.path.join(outdir, "input_data_used.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dataset", "species", "category", "category_tulv_moved",
                    "host_class", "is_rodent", "NW_HPS_clade", "n_sequences",
                    "piN", "piS", "piN_piS", "flag_piS_saturated"])
        for ds in out["datasets"]:
            if ds.get("skipped"):
                continue
            for s_ in ds["species"]:
                w.writerow([ds["tag"], s_, ds["category"][s_],
                            "pathogenic" if (ds["category"][s_] == "environmental"
                                             and s_ == "TULV") else ds["category"][s_],
                            ds["host_class"][s_], ds["is_rodent"][s_],
                            ds["NW_HPS_clade"][s_], ds["n_sequences_used"][s_],
                            ds["piN"][s_], ds["piS"][s_], ds["piN_piS"][s_],
                            ds["saturation_flags"][s_]])
    print("wrote", os.path.join(outdir, "input_data_used.tsv"))
    print("wrote", os.path.join(outdir, "phylo_stats.json"))
    print("wrote", os.path.join(outdir, "pgls_model_fits.tsv"))


if __name__ == "__main__":
    main()
