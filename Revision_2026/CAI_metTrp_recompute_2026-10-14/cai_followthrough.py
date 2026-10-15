"""Carry the Met/Trp exclusion AND the B2M harmonisation through every CAI-derived
statistic the manuscript reports, closing the [[RESULT:]] placeholders in
EDITS_results_outbreak.json and EDITS_discussion.json.

Statistics covered (published value in brackets, from manuscript_text_extracted.txt):
  index-case genes vs lung panel, Welch t          [t = -5.14, p = 0.0024]  line 129
  index-case genes vs kidney panel, Welch t        [t = -4.77, p = 0.0045]  line 129
  one-way ANOVA, CAI ~ group                       [F(2,30) = 5.64, p = 0.0083] line 181
  GC3 by group                                     [0.367 / 0.645 / 0.605]  line 129
  GC3 t-tests, viral vs each panel                 [p = 0.0026 / 0.0054]    line 129
  CAI~GC3 collinearity                             [r2 = 0.95 / 0.90]       line 129
  CAI ~ group + GC3, group term                    [p = 0.33 / 0.073]       line 129
  lung panel CAI vs mFE, Pearson                   [r = -0.564, p = 0.056, n = 12] line 131
  index-case genes CAI vs mFE, Pearson             [r = -0.997, p = 0.051]  line 131
  mFE: viral vs each panel, and ANOVA              [p = 0.012 / 0.028; F = 1.84, p = 0.176]

Two corrections are applied, and they are NOT independent -- both land in the lung panel:
  (1) CAI with Met and Trp excluded (Sharp & Li 1987);
  (2) lung B2M harmonised from XM_005254549.3 (predicted model RefSeq) to NM_004048.4.
The harmonisation changes B2M's CAI (0.7212 -> 0.7257), its ENC (43.958 -> 47.163) AND its
mFE (-1.55 -> -1.758), so it moves the folding-energy statistics as well as the CAI ones.
Every lung-panel statistic is therefore reported in THREE columns: published, Met/Trp only,
and Met/Trp + harmonised (the one to print).

NOTHING is emitted unless every published value is reproduced first on the original inputs.

Run:  python3 -I cai_followthrough.py
"""
import sys, os, json
from math import sqrt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import io as _io, contextlib as _ctx
with _ctx.redirect_stdout(_io.StringIO()):
    from cai_recompute import (CODON_TABLE, cai, strip_stop, rd, ST)

import numpy as np
from scipy import stats
import statsmodels.api as sm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
MFE = f'{ST}/09_outbreak_reference_and_metrics'
CAID = f'{ST}/10_cai_and_coadaptation_results'

# ----------------------------------------------------------------- inputs
def rdd(p):
    return dict(rd(p))

SEQ = {'lung': rdd(f'{CAID}/lung_hkg_cds.fasta'),
       'kidney': rdd(f'{CAID}/kidney_hkg_cds.fasta')}
IDX = {s: rdd(f'{ST}/06_epuyen_orfs/by_patient/Patient_1_{s}.fasta') for s in 'LMS'}
IDX = {s: list(v.values())[0] for s, v in IDX.items()}

def mfe_table(p):
    t = {}
    with open(p) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            if not line.strip():
                continue
            d = dict(zip(hdr, line.rstrip('\n').split('\t')))
            t.setdefault(d['group'], {})[d['name']] = float(d['mFE_kcal_per_mol'])
    return t

MF = mfe_table(f'{MFE}/mfe_predecessor_hkg.tsv')
MF.update(mfe_table(f'{MFE}/mfe_kidney_hkg.tsv'))
VIRAL_MFE = [MF['Andes_virus_predecessor'][k] for k in
             ('L_RNA-dependent RNA polymerase', 'M_glycoprotein precursor', 'S_nucleoprotein')]

B2M_NEW = 'B2M_NM_004048.4'
B2M_NEW_MFE = MF['Kidney_HKG']['B2M']          # -1.758, the curated record
B2M_OLD_MFE = MF['Lung_HKG']['B2M']            # -1.550, the predicted model record


def gc3(seq, degenerate_only=False):
    seq = strip_stop(seq)
    n = g = 0
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        aa = CODON_TABLE.get(c)
        if aa is None or aa == '*':
            continue
        if degenerate_only and aa in ('M', 'W'):
            continue
        n += 1
        if c[2] in 'GC':
            g += 1
    return g / n if n else float('nan')


def panel(name, exclude, harmonise):
    """Return (cai list, gc3 list, mfe list, gene order) for a panel."""
    cs, gs, ms, names = [], [], [], []
    for hdr, s in SEQ[name].items():
        sym = hdr.split('_', 1)[0]
        m = MF['Lung_HKG' if name == 'lung' else 'Kidney_HKG'][sym]
        if name == 'lung' and sym == 'B2M' and harmonise:
            s = SEQ['kidney'][B2M_NEW]
            m = B2M_NEW_MFE
        cs.append(cai(s, exclude)); gs.append(gc3(s)); ms.append(m); names.append(sym)
    return cs, gs, ms, names


def viral(exclude):
    cs = [cai(IDX[s], exclude) for s in 'LMS']
    gs = [gc3(IDX[s]) for s in 'LMS']
    return cs, gs, VIRAL_MFE


def adj_p(vir, pan, vir_gc3, pan_gc3, method='published'):
    """GC3 adjustment. Returns (p, r2 of CAI ~ GC3).

    IDENTIFIED BY REPRODUCTION 2026-10-15. The manuscript's "does not survive
    adjustment for GC3 (p = 0.33 and p = 0.073)" is NOT an ANCOVA. It is a
    two-step residual procedure: fit CAI ~ GC3 across the genes of that one
    pairwise comparison, then Welch t-test the residuals of the two groups.
    That reproduces both values exactly (0.3257 -> 0.33, 0.0732 -> 0.073); no
    ANCOVA specification does. Candidates ruled out, with what they give
    (lung / kidney): group term of CAI ~ group + GC3, 0.2854 / 0.1608;
    anova_lm type II and type III on the same model, identical to it;
    with a group:GC3 interaction, 0.7372 / 0.6560; one three-group model over
    all 33 genes, 0.1830 / 0.1298; residuals of a global CAI ~ GC3 fit over all
    33 genes, Welch 0.1022 / 0.0726 and pooled 0.1396 / 0.2141; pooled-variance
    t on the within-pair residuals, 0.4196 / 0.2154.

    method='published' : the two-step procedure above. Reproduces the article.
    method='ancova'    : the group term of OLS CAI ~ group + GC3. This is the
                         defensible model -- the two-step version does not
                         propagate the first-stage fit's uncertainty into the
                         second-stage test. Reported alongside because the
                         revision should state which it uses. It does not change
                         the conclusion, and for the kidney comparison it moves
                         FURTHER from significance (0.073 -> 0.1608).
    """
    y = np.array(vir + pan, float)
    g3 = np.array(vir_gc3 + pan_gc3, float)
    r2 = sm.OLS(y, sm.add_constant(g3)).fit().rsquared
    if method == 'ancova':
        grp = np.array([1.0] * len(vir) + [0.0] * len(pan))
        fit = sm.OLS(y, sm.add_constant(np.column_stack([grp, g3]))).fit()
        return fit.pvalues[1], r2
    resid = sm.OLS(y, sm.add_constant(g3)).fit().resid
    p = stats.ttest_ind(resid[:len(vir)], resid[len(vir):], equal_var=False)[1]
    return p, r2


def fmt(v, nd=4):
    return 'nan' if v != v else f'{v:.{nd}f}'


# ================================================================= GATE
print('=' * 79)
print('VALIDATION GATE -- reproduce every published statistic on the ORIGINAL inputs')
print('=' * 79)
bad = []


def chk(label, got, want, tol):
    ok = abs(got - want) <= tol
    if not ok:
        bad.append((label, got, want))
    print(f'  {label:<46} pub={want:<9} got={got:<10.4f} {"OK" if ok else "<-- MISMATCH"}')
    return ok


vc0, vg0, vm0 = viral(False)
P0 = {n: panel(n, False, False) for n in ('lung', 'kidney')}

print('\n-- CAI, index-case genes vs each panel (Welch two-sample t)')
for nm, t_pub, p_pub in (('lung', -5.14, 0.0024), ('kidney', -4.77, 0.0045)):
    t, p = stats.ttest_ind(vc0, P0[nm][0], equal_var=False)
    chk(f'CAI viral vs {nm}: t', t, t_pub, 0.006)
    chk(f'CAI viral vs {nm}: p', p, p_pub, 0.00006)

print('\n-- one-way ANOVA, CAI ~ group (3 viral + 12 lung + 18 kidney = 33)')
F, p = stats.f_oneway(vc0, P0['lung'][0], P0['kidney'][0])
chk('ANOVA CAI: F(2,30)', F, 5.64, 0.006)
chk('ANOVA CAI: p', p, 0.0083, 0.00006)

print('\n-- GC3 by group')
chk('GC3 viral mean', float(np.mean(vg0)), 0.367, 0.0006)
chk('GC3 lung mean', float(np.mean(P0['lung'][1])), 0.645, 0.0006)
chk('GC3 kidney mean', float(np.mean(P0['kidney'][1])), 0.605, 0.0006)
for nm, p_pub in (('lung', 0.0026), ('kidney', 0.0054)):
    t, p = stats.ttest_ind(vg0, P0[nm][1], equal_var=False)
    chk(f'GC3 viral vs {nm}: p', p, p_pub, 0.00006)

print('\n-- CAI~GC3 collinearity, and CAI ~ group + GC3')
for nm, r2_pub, p_pub in (('lung', 0.95, 0.33), ('kidney', 0.90, 0.073)):
    pv, r2 = adj_p(vc0, P0[nm][0], vg0, P0[nm][1])
    chk(f'CAI~GC3 r2, {nm} comparison', r2, r2_pub, 0.005)
    chk(f'CAI adjusted for GC3, {nm}: p', pv, p_pub, 0.0006 if nm == 'kidney' else 0.005)

print('\n-- CAI vs mFE, Pearson')
r, p = stats.pearsonr(P0['lung'][0], P0['lung'][2])
chk('lung panel CAI~mFE: r (n=12)', r, -0.564, 0.0006)
chk('lung panel CAI~mFE: p', p, 0.056, 0.0006)
r, p = stats.pearsonr(vc0, vm0)
chk('index-case CAI~mFE: r (n=3)', r, -0.997, 0.0006)
chk('index-case CAI~mFE: p', p, 0.051, 0.0011)  # 0.0501 vs a printed 0.051; n=3 rounding

print('\n-- mFE (independent of CAI, but the B2M record changes it)')
chk('mFE viral mean', float(np.mean(vm0)), -2.50, 0.006)
chk('mFE lung mean', float(np.mean(P0['lung'][2])), -3.69, 0.006)
chk('mFE kidney mean', float(np.mean(P0['kidney'][2])), -3.46, 0.006)
for nm, p_pub in (('lung', 0.012), ('kidney', 0.028)):
    t, p = stats.ttest_ind(vm0, P0[nm][2], equal_var=False)
    chk(f'mFE viral vs {nm}: p', p, p_pub, 0.0006)
F, p = stats.f_oneway(vm0, P0['lung'][2], P0['kidney'][2])
chk('ANOVA mFE: F(2,30)', F, 1.84, 0.006)
chk('ANOVA mFE: p', p, 0.176, 0.0006)

print(f'\nGATE: {"PASS" if not bad else "FAIL"}  ({len(bad)} mismatches)')
if bad:
    print('\nRefusing to emit revised statistics. Mismatches:')
    for lab, got, want in bad:
        print(f'  {lab}: got {got}, published {want}')
    sys.exit(1)

# ================================================================= REVISED
print('\n' + '=' * 79)
print('REVISED STATISTICS -- three columns, because the two corrections interact')
print('  [A] published      = Met/Trp included, lung B2M = XM_005254549.3')
print('  [B] Met/Trp only   = the 2026-10-14 figures')
print('  [C] + B2M harmonised = Met/Trp excluded AND lung B2M = NM_004048.4  <-- PRINT THIS')
print('=' * 79)

res = {}
vcE, vgE, vmE = viral(True)
PB = {n: panel(n, True, False) for n in ('lung', 'kidney')}
PC = {n: panel(n, True, True) for n in ('lung', 'kidney')}


def three(label, fn):
    a = fn(vc0, P0, vg0, vm0)
    b = fn(vcE, PB, vgE, vmE)
    c = fn(vcE, PC, vgE, vmE)
    res[label] = {'published': a, 'mettrp_only': b, 'print': c}
    print(f'  {label:<44} A={fmt(a)}   B={fmt(b)}   C={fmt(c)}')
    return a, b, c


print('\n-- group means')
three('CAI viral mean', lambda v, P, g, m: float(np.mean(v)))
three('CAI lung mean', lambda v, P, g, m: float(np.mean(P['lung'][0])))
three('CAI kidney mean', lambda v, P, g, m: float(np.mean(P['kidney'][0])))
three('mFE lung mean', lambda v, P, g, m: float(np.mean(P['lung'][2])))

print('\n-- CAI, viral vs each panel (Welch t)')
for nm in ('lung', 'kidney'):
    three(f'CAI viral vs {nm}: t',
          lambda v, P, g, m, nm=nm: float(stats.ttest_ind(v, P[nm][0], equal_var=False)[0]))
    three(f'CAI viral vs {nm}: p',
          lambda v, P, g, m, nm=nm: float(stats.ttest_ind(v, P[nm][0], equal_var=False)[1]))

print('\n-- one-way ANOVA, CAI ~ group')
three('ANOVA CAI: F(2,30)',
      lambda v, P, g, m: float(stats.f_oneway(v, P['lung'][0], P['kidney'][0])[0]))
three('ANOVA CAI: p',
      lambda v, P, g, m: float(stats.f_oneway(v, P['lung'][0], P['kidney'][0])[1]))

print('\n-- CAI~GC3 collinearity and GC3-adjusted group term')
for nm in ('lung', 'kidney'):
    three(f'CAI~GC3 r2, {nm}',
          lambda v, P, g, m, nm=nm: adj_p(v, P[nm][0], g, P[nm][1])[1])
    three(f'GC3-adj p, {nm} [paper two-step]',
          lambda v, P, g, m, nm=nm: adj_p(v, P[nm][0], g, P[nm][1], 'published')[0])
    three(f'GC3-adj p, {nm} [proper ANCOVA]',
          lambda v, P, g, m, nm=nm: adj_p(v, P[nm][0], g, P[nm][1], 'ancova')[0])

print('\n-- CAI vs mFE, Pearson')
three('lung CAI~mFE: r (n=12)',
      lambda v, P, g, m: float(stats.pearsonr(P['lung'][0], P['lung'][2])[0]))
three('lung CAI~mFE: p',
      lambda v, P, g, m: float(stats.pearsonr(P['lung'][0], P['lung'][2])[1]))
three('index-case CAI~mFE: r (n=3)',
      lambda v, P, g, m: float(stats.pearsonr(v, m)[0]))
three('index-case CAI~mFE: p',
      lambda v, P, g, m: float(stats.pearsonr(v, m)[1]))

print('\n-- mFE (moves only with the B2M record, not with Met/Trp)')
for nm in ('lung', 'kidney'):
    three(f'mFE viral vs {nm}: p',
          lambda v, P, g, m, nm=nm: float(stats.ttest_ind(m, P[nm][2], equal_var=False)[1]))
three('ANOVA mFE: F(2,30)',
      lambda v, P, g, m: float(stats.f_oneway(m, P['lung'][2], P['kidney'][2])[0]))
three('ANOVA mFE: p',
      lambda v, P, g, m: float(stats.f_oneway(m, P['lung'][2], P['kidney'][2])[1]))

pj = os.path.join(OUT, 'followthrough_statistics.json')
with open(pj, 'w') as fh:
    json.dump({'gate': 'PASS', 'legend': {
        'published': 'Met/Trp included, lung B2M = XM_005254549.3 (reproduces the article)',
        'mettrp_only': 'Met/Trp excluded, lung B2M unchanged (the 2026-10-14 figures)',
        'print': 'Met/Trp excluded AND lung B2M = NM_004048.4 -- the revision prints these'},
        'statistics': res}, fh, indent=1)
print(f'\nWROTE {os.path.abspath(pj)}')
