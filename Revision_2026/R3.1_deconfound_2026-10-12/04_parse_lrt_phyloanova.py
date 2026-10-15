#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""R3.1 — step 4: collect every codeml fit, run the nested LRTs, and run a
phylogenetic ANOVA on the per-species omega estimates.

Branch models, all on the SAME tree and SAME 448-codon alignment, terminal
branches labelled only (ancestral host/disease states unknown):

  M0        one omega for the whole tree                       (from R2.2)
  disease2  foreground = 10 pathogenic tips                    (from R2.2)
  host2     foreground = 12 rodent-borne tips                  (here)
  three     w1 = pathogenic tips, w2 = PHV+TULV, w0 = rest     (here)
  cladeNW   foreground = New World HPS clade + its internals   (here)
  disease2_tulvpath  R3.3c: TULV counted pathogenic            (here)
  disease2_tulvown   R3.3c: TULV given its own omega           (here)

Nesting:  disease2 = three with w2 == w0 ;  host2 = three with w1 == w2.
So  three vs host2   (df 1) = DISEASE effect *within* the rodent-borne host
                              class  -- this is R3.1a done without discarding
                              the rest of the tree;
    three vs disease2 (df 1) = HOST effect *within* the environmental class.
disease2 and host2 have identical np, so their lnL compare directly (R3.1b).
"""
import glob
import json
import os
import re
import sys

import numpy as np
from scipy import stats
from Bio import Phylo

A = ("/PATH/TO/"
     "Project_Hantavirus/Reviews_ans/analysis/R3.1_deconfound_2026-10-12")
R22 = ("/PATH/TO/"
       "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12")

PATHOGENIC = ["HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV",
              "CHOV", "LANV"]
RODENT = PATHOGENIC + ["PHV", "TULV"]
ALL22 = RODENT + ["TPMV", "MJNV", "NVAV", "ASAV", "OXBV", "RKPV", "SWSV",
                  "CBNV", "QZNV", "JJUV"]
HOSTCLASS = {}
for t in ALL22:
    HOSTCLASS[t] = "rodent" if t in RODENT else "non-rodent"


def parse(path):
    """Pull lnL, np, kappa and the branch omegas out of a codeml out.txt."""
    if not os.path.exists(path):
        return None
    txt = open(path).read()
    m = re.search(r"^lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\):\s*(-?[\d.]+)",
                  txt, re.M)
    if not m:
        return None
    d = dict(ntime=int(m.group(1)), np=int(m.group(2)), lnL=float(m.group(3)))
    k = re.search(r"^kappa \(ts/tv\) =\s*([\d.]+)", txt, re.M)
    d["kappa"] = float(k.group(1)) if k else None
    w = re.search(r"^w \(dN/dS\) for branches:\s*(.+)$", txt, re.M)
    if w:
        d["w"] = [float(x) for x in w.group(1).split()]
    else:
        o = re.search(r"^omega \(dN/dS\) =\s*([\d.]+)", txt, re.M)
        d["w"] = [float(o.group(1))] if o else None
    return d


fits = {}
for cf in ("F3x4", "F61"):
    fits[f"M0_{cf}"] = parse(f"{R22}/results/codeml/m0_{cf.lower() if cf=='F61' else 'F3x4'}/out.txt")
# the R2.2 directory names are m0_F3x4 / m0_F61 / two_ratio_F3x4 / two_ratio_F61
fits["M0_F3x4"] = parse(f"{R22}/results/codeml/m0_F3x4/out.txt")
fits["M0_F61"] = parse(f"{R22}/results/codeml/m0_F61/out.txt")
fits["disease2_F3x4"] = parse(f"{R22}/results/codeml/two_ratio_F3x4/out.txt")
fits["disease2_F61"] = parse(f"{R22}/results/codeml/two_ratio_F61/out.txt")
for d in sorted(glob.glob(f"{A}/results/codeml/*/out.txt")):
    name = os.path.basename(os.path.dirname(d))
    fits[name] = parse(d)

out = []
def say(s=""):
    out.append(s)
    print(s)

say("# R3.1 — branch-model head-to-head: disease labelling vs host labelling")
say()
say("Same ML tree (IQ-TREE 3.1.1, GTR+F+R4), same 22-species / 448-codon S-segment")
say("alignment, identical codeml control files, terminal branches labelled only.")
say("`M0` and `disease2` are the runs already in")
say("`R2.2_phylo_dnds_2026-10-12/results/codeml/` (m0_* and two_ratio_*);")
say("everything else was fitted here.")
say()

say("## Fitted models")
say()
say("| model | codon freq | np | lnL | kappa | omegas (background, then #1, #2) |")
say("|---|---|---|---|---|---|")
order = ["M0", "disease2", "host2", "three", "cladeNW",
         "disease2_tulvpath", "disease2_tulvown"]
for cf in ("F3x4", "F61"):
    for base in order:
        key = f"{base}_{cf}"
        f = fits.get(key)
        if f is None:
            say(f"| {base} | {cf} | — | **PENDING / not found** | — | — |")
            continue
        ws = " ".join(f"{x:.5f}" for x in f["w"])
        say(f"| {base} | {cf} | {f['np']} | {f['lnL']:.6f} | "
            f"{f['kappa']:.5f} | {ws} |")
say()

say("## Which codon-frequency model should be believed?")
say()
say("This matters because the two codon-frequency models disagree about the")
say("SIGN of the pathogenic-vs-background omega contrast, so the choice is not")
say("cosmetic. codeml's reported `np` excludes the codon-frequency parameters,")
say("which PAML estimates by counting rather than by ML: F3x4 (CodonFreq=2)")
say("uses 9 free nucleotide frequencies, F61 (CodonFreq=3) uses 60 free codon")
say("frequencies. Counting those gives the honest penalty.")
say()
say("| M0 fit | lnL | np (codeml) | + frequency params | total | AIC |")
say("|---|---|---|---|---|---|")
_freq = {"F3x4": 9, "F61": 60}
_aic = {}
for cf in ("F3x4", "F61"):
    f = fits.get(f"M0_{cf}")
    if not f:
        continue
    tot = f["np"] + _freq[cf]
    aic = 2 * (-f["lnL"]) + 2 * tot
    _aic[cf] = aic
    say(f"| M0 {cf} | {f['lnL']:.4f} | {f['np']} | {_freq[cf]} | {tot} | "
        f"{aic:.2f} |")
say()
if len(_aic) == 2:
    f3, f6 = fits["M0_F3x4"], fits["M0_F61"]
    d = 2 * (f6["lnL"] - f3["lnL"])
    ddf = _freq["F61"] - _freq["F3x4"]
    say(f"2dlnL(F61 vs F3x4) = {d:.2f} on df = {ddf}, "
        f"p = {stats.chi2.sf(d, ddf):.3g}; "
        f"dAIC = {_aic['F3x4'] - _aic['F61']:.2f} in favour of **F61**.")
    say()
    say("**F61 is decisively the better-fitting codon-frequency model**, and it")
    say("is also the model this comparison calls for (\"codeml M0 with a")
    say("codon-frequency model\"). Any conclusion that holds only under F3x4 is")
    say("therefore the weaker of the two, and must be reported that way.")
    say()


def lrt(alt, null, df, label):
    a, n = fits.get(alt), fits.get(null)
    if a is None or n is None:
        return f"| {label} | — | — | — | **PENDING** |"
    d = 2 * (a["lnL"] - n["lnL"])
    p = stats.chi2.sf(d, df) if d > 0 else 1.0
    return (f"| {label} | {a['lnL']:.4f} vs {n['lnL']:.4f} | {d:.3f} | {df} | "
            f"{p:.4g} |")

say("## Likelihood-ratio tests")
say()
for cf in ("F3x4", "F61"):
    say(f"### {cf}")
    say()
    say("| test | lnL alt vs null | 2dlnL | df | p |")
    say("|---|---|---|---|---|")
    say(lrt(f"disease2_{cf}", f"M0_{cf}", 1,
            "disease labelling vs one omega"))
    say(lrt(f"host2_{cf}", f"M0_{cf}", 1,
            "host labelling vs one omega"))
    say(lrt(f"cladeNW_{cf}", f"M0_{cf}", 1,
            "New World HPS clade vs one omega"))
    say(lrt(f"three_{cf}", f"M0_{cf}", 2,
            "three classes vs one omega"))
    say(lrt(f"three_{cf}", f"host2_{cf}", 1,
            "**DISEASE effect within rodent-borne** (three vs host)"))
    say(lrt(f"three_{cf}", f"disease2_{cf}", 1,
            "**HOST effect within environmental** (three vs disease)"))
    say(lrt(f"disease2_tulvown_{cf}", f"disease2_{cf}", 1,
            "R3.3c: TULV given its own omega vs TULV as environmental"))
    say()
    # non-nested, equal np: direct lnL / AIC comparison
    d, h = fits.get(f"disease2_{cf}"), fits.get(f"host2_{cf}")
    c = fits.get(f"cladeNW_{cf}")
    tp = fits.get(f"disease2_tulvpath_{cf}")
    if d and h:
        say(f"**Head-to-head (equal np = {d['np']}, so lnL compares directly):**")
        say()
        say("| labelling of the SAME tree | lnL | dAIC vs best |")
        say("|---|---|---|")
        cand = {"disease (10 pathogenic tips)": d,
                "host (12 rodent-borne tips)": h}
        if c:
            cand["clade (New World HPS clade)"] = c
        if tp:
            cand["disease, TULV counted pathogenic"] = tp
        best = max(x["lnL"] for x in cand.values())
        for nm, f in sorted(cand.items(), key=lambda kv: -kv[1]["lnL"]):
            say(f"| {nm} | {f['lnL']:.4f} | {2 * (best - f['lnL']):.3f} |")
        say()
        win = max(cand.items(), key=lambda kv: kv[1]["lnL"])[0]
        say(f"Best-fitting labelling under {cf}: **{win}**.")
        say()

# ------------------------------------------------------------------ R3.1a ----
say("## R3.1a — the DECONFOUNDED disease contrast: pathogenic vs environmental")
say("## *within* the rodent-borne host class")
say()
say("In the three-class model `w1` is the 10 pathogenic terminal branches and")
say("`w2` is the 2 environmental rodent-borne terminal branches (PHV, TULV);")
say("host class is identical for both, so this contrast is not confounded by")
say("host. It is, however, 10 species against 2.")
say()
say("| codon freq | w1 pathogenic | w2 rodent-environmental | ratio w2/w1 | LRT vs w1==w2 | p | direction |")
say("|---|---|---|---|---|---|---|")
for cf in ("F3x4", "F61"):
    t, h = fits.get(f"three_{cf}"), fits.get(f"host2_{cf}")
    if not (t and h and t.get("w") and len(t["w"]) >= 3):
        say(f"| {cf} | — | — | — | **PENDING** | | |")
        continue
    w1, w2 = t["w"][1], t["w"][2]
    d = 2 * (t["lnL"] - h["lnL"])
    p = stats.chi2.sf(d, 1) if d > 0 else 1.0
    direc = ("pathogenic LOWER (supports the manuscript)" if w1 < w2
             else "pathogenic HIGHER (contradicts the manuscript)")
    say(f"| {cf} | {w1:.5f} | {w2:.5f} | {w2 / w1:.2f} | {d:.3f} | {p:.4g} | "
        f"{direc} |")
say()
say("### Why the two-class disease model appeared to reverse under F61")
say()
say("The two-class disease model forces PHV and TULV into the same background")
say("class as the ten non-rodent environmental species. The three-class fits")
say("show those two groups have very different omega, so that background class")
say("is a mixture and the two-class contrast is misspecified:")
say()
say("| codon freq | w0 background (non-rodent env + internal branches) | w1 pathogenic tips | w2 PHV+TULV |")
say("|---|---|---|---|")
for cf in ("F3x4", "F61"):
    t = fits.get(f"three_{cf}")
    if t and t.get("w") and len(t["w"]) >= 3:
        say(f"| {cf} | {t['w'][0]:.5f} | {t['w'][1]:.5f} | {t['w'][2]:.5f} |")
say()

say("### Does the deconfounded contrast rest on a single species?")
say()
say("PHV and TULV are a sister pair on this tree (`{PHV,TULV}` is a clade, see")
say("`01_topology_check.md`), so as a contrast they amount to ONE")
say("phylogenetically independent observation. Repeating the within-rodent test")
say("with each of them alone:")
say()
say("| codon freq | environmental rodent species used | its omega | pathogenic omega | 2dlnL | p |")
say("|---|---|---|---|---|---|")
for cf in ("F3x4", "F61"):
    for who, altk, nullk in (
            ("PHV only", f"three_phvonly_{cf}", f"disease2_phvpath_{cf}"),
            ("TULV only", f"disease2_tulvown_{cf}", f"disease2_tulvpath_{cf}"),
            ("PHV + TULV", f"three_{cf}", f"host2_{cf}")):
        a, n = fits.get(altk), fits.get(nullk)
        if not (a and n and a.get("w") and len(a["w"]) >= 3):
            say(f"| {cf} | {who} | — | — | **PENDING** | |")
            continue
        d = 2 * (a["lnL"] - n["lnL"])
        p = stats.chi2.sf(d, 1) if d > 0 else 1.0
        say(f"| {cf} | {who} | {a['w'][2]:.5f} | {a['w'][1]:.5f} | {d:.3f} | "
            f"{p:.4g} |")
say()

# ------------------------------------------------------------------ R3.1c ----
say("## R3.1c — per-species omega and a phylogenetically controlled comparison")
say()
say("Each species in turn was made the sole foreground branch (22 codeml runs")
say("per codon-frequency model), giving a per-species omega estimated against")
say("the whole-tree background. Boundary hits (omega at 0 or at codeml's upper")
say("limit) are flagged, because they are not real estimates.")
say()

tree = Phylo.read(f"{A}/inputs/S22_ml_clean.nwk", "newick")
names = [t.name for t in tree.get_terminals()]
# variance-covariance matrix under Brownian motion = root-to-MRCA distance.
# The IQ-TREE tree has no outgroup, so the root placement is arbitrary; this is
# stated as a caveat rather than hidden.
root = tree.root
vcv = np.zeros((22, 22))
idx = {n: i for i, n in enumerate(names)}
for i, a in enumerate(names):
    for j, b in enumerate(names):
        if i == j:
            vcv[i, j] = tree.distance(root, tree.find_any(a))
        else:
            mrca = tree.common_ancestor(a, b)
            vcv[i, j] = tree.distance(root, mrca)

def phyloanova(trait, groups, nsim=20000, seed=20261012):
    """Garland et al. (1993) phylogenetic ANOVA: observed F vs an F null built
    from Brownian-motion simulations on the tree."""
    trait = np.asarray(trait, float)
    g = np.asarray(groups)
    levels = sorted(set(g))
    def F(y):
        gm = y.mean()
        ssb = sum(len(y[g == L]) * (y[g == L].mean() - gm) ** 2 for L in levels)
        ssw = sum(((y[g == L] - y[g == L].mean()) ** 2).sum() for L in levels)
        dfb, dfw = len(levels) - 1, len(y) - len(levels)
        if ssw <= 0 or dfw <= 0:
            return np.nan
        return (ssb / dfb) / (ssw / dfw)
    f_obs = F(trait)
    rng = np.random.default_rng(seed)
    # Cholesky of the VCV; add a tiny ridge for numerical positive-definiteness
    L = np.linalg.cholesky(vcv + 1e-9 * np.eye(22))
    sims = L @ rng.standard_normal((22, nsim))
    sims *= trait.std(ddof=1) / sims.std(ddof=1)
    fs = np.array([F(sims[:, k]) for k in range(nsim)])
    fs = fs[~np.isnan(fs)]
    p_phy = float((np.sum(fs >= f_obs) + 1) / (len(fs) + 1))
    dfb, dfw = len(levels) - 1, len(trait) - len(levels)
    p_std = float(stats.f.sf(f_obs, dfb, dfw))
    return dict(F=float(f_obs), p_standard=p_std, p_phylo=p_phy, nsim=len(fs))

pertip = {}
for cf in ("F3x4", "F61"):
    row = {}
    for sp in ALL22:
        f = fits.get(f"pertip_{sp}_{cf}")
        row[sp] = f["w"][1] if (f and f["w"] and len(f["w"]) > 1) else None
    pertip[cf] = row

for cf in ("F3x4", "F61"):
    row = pertip[cf]
    miss = [s for s, v in row.items() if v is None]
    if miss:
        say(f"### {cf} — **PENDING**: {len(miss)} of 22 per-species runs not "
            f"finished ({', '.join(miss)}). Values below omit them.")
        say()
    got = {s: v for s, v in row.items() if v is not None}
    if len(got) < 12:
        say(f"Too few {cf} estimates to test; rerun script 04 when complete.")
        say()
        continue
    say(f"### {cf} — per-species omega (foreground = that species' own branch)")
    say()
    say("| species | host | category | omega | boundary? |")
    say("|---|---|---|---|---|")
    for sp in ALL22:
        v = row[sp]
        if v is None:
            say(f"| {sp} | {HOSTCLASS[sp]} | "
                f"{'pathogenic' if sp in PATHOGENIC else 'environmental'} | "
                f"PENDING | |")
            continue
        b = "yes (lower bound)" if v <= 1e-5 else (
            "yes (upper bound)" if v >= 998 else "")
        say(f"| {sp} | {HOSTCLASS[sp]} | "
            f"{'pathogenic' if sp in PATHOGENIC else 'environmental'} | "
            f"{v:.5f} | {b} |")
    say()
    sp_ok = [s for s in ALL22 if row[s] is not None]
    nb = [s for s in sp_ok if 1e-5 < row[s] < 998]
    say(f"Usable (non-boundary) estimates: {len(nb)} of {len(sp_ok)}.")
    say()
    if len(nb) >= 12:
        y = np.array([row[s] for s in nb])
        ylog = np.log(y)
        full = [s for s in names]           # tree order for the VCV
        # restrict the VCV to the usable species, keeping tree order
        keep = [i for i, n in enumerate(names) if n in nb]
        sub = vcv[np.ix_(keep, keep)]
        trait = np.array([row[names[i]] for i in keep])
        traitlog = np.log(trait)
        gdis = np.array(["pathogenic" if names[i] in PATHOGENIC
                         else "environmental" for i in keep])
        ghos = np.array([HOSTCLASS[names[i]] for i in keep])
        # temporarily swap in the restricted VCV
        vcv_full = vcv
        globals()["vcv"] = sub
        say("| trait | grouping | mean group 1 | mean group 2 | F | p (standard ANOVA) | p (phylogenetic ANOVA) |")
        say("|---|---|---|---|---|---|---|")
        for tname, tv in (("omega", trait), ("log omega", traitlog)):
            for gname, gv in (("disease", gdis), ("host", ghos)):
                r = phyloanova(tv, gv)
                L = sorted(set(gv))
                m1 = tv[gv == L[0]].mean()
                m2 = tv[gv == L[1]].mean()
                say(f"| {tname} | {gname} | {L[0]} {m1:.4f} | {L[1]} {m2:.4f} | "
                    f"{r['F']:.3f} | {r['p_standard']:.4g} | "
                    f"**{r['p_phylo']:.4g}** |")
        globals()["vcv"] = vcv_full
        say()
        say("The phylogenetic-ANOVA p is the fraction of 20 000 Brownian-motion")
        say("simulations on this tree that produce an F at least as large as the")
        say("observed one. It is the correct null for 'could the tree alone")
        say("produce this group difference?'. Caveat: the tree has no outgroup,")
        say("so the root placement used for the variance-covariance matrix is")
        say("the arbitrary IQ-TREE rooting.")
        say()

with open(f"{A}/results/04_branch_models_and_phyloanova.md", "w") as fh:
    fh.write("\n".join(out) + "\n")
with open(f"{A}/results/04_fits.json", "w") as fh:
    json.dump({"fits": fits, "pertip": pertip}, fh, indent=1)
print("\nwrote results/04_branch_models_and_phyloanova.{md,json}")
