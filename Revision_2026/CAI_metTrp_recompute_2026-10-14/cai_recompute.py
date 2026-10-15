# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""Recompute CAI with Met/Trp EXCLUDED (Sharp & Li 1987) alongside the
published INCLUDED-at-w=1 behaviour. Include-mode must reproduce the paper."""
import sys, glob, os
from math import log, exp
from statistics import mean

B = '/PATH/TO/Project_Hantavirus'
ST = B + '/Reviews_ans/zenodo/staging'

BASES = 'TCAG'
AAS = ('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
CODON_TABLE = {}
i = 0
for a in BASES:
    for b in BASES:
        for c in BASES:
            CODON_TABLE[a+b+c] = AAS[i]; i += 1

AA_TO_CODONS = {}
for c, aa in CODON_TABLE.items():
    if aa != '*':
        AA_TO_CODONS.setdefault(aa, []).append(c)

SINGLE = {c for aa, cs in AA_TO_CODONS.items() if len(cs) == 1 for c in cs}

def weights():
    raw = {}
    for line in open(B + '/Work1/Epuyen_outbreak/Ehum.cut'):
        line = line.strip()
        if not line or line.startswith('#'): continue
        p = line.split()
        if p[1] == '*': continue
        raw[p[0].upper().replace('U','T')] = float(p[4])
    w = {}
    for aa, cs in AA_TO_CODONS.items():
        fam = {c: raw.get(c, 0.0) for c in cs}
        mx = max(fam.values())
        for c, f in fam.items():
            w[c] = (f/mx) if mx > 0 else 1.0
    return w

W = weights()

def strip_stop(s):
    s = s.upper().replace('U','T')
    if len(s) % 3 == 0 and s[-3:] in ('TAA','TAG','TGA'): return s[:-3]
    return s

def cai(seq, exclude_single):
    seq = strip_stop(seq)
    ls, n = 0.0, 0
    for i in range(0, len(seq)-2, 3):
        c = seq[i:i+3]
        if c not in W: continue
        if exclude_single and c in SINGLE: continue
        ls += log(max(W[c], 1e-6)); n += 1
    return exp(ls/n) if n else float('nan')

def rd(p):
    recs, h, s = [], None, []
    for line in open(p):
        line = line.rstrip()
        if not line: continue
        if line.startswith('>'):
            if h is not None: recs.append((h, ''.join(s)))
            h, s = line[1:], []
        else: s.append(line)
    if h is not None: recs.append((h, ''.join(s)))
    return recs

def row(label, vals_inc, vals_exc, published=None):
    mi, me = mean(vals_inc), mean(vals_exc)
    pub = f'{published:.4f}' if published is not None else '  --  '
    ok = ''
    if published is not None:
        ok = ' MATCH' if abs(mi - published) < 0.0006 else f' <-- MISMATCH ({mi-published:+.4f})'
    print(f'  {label:<34} pub={pub}  inc={mi:.4f}  exc={me:.4f}  d={me-mi:+.4f}{ok}')
    return mi, me

print('Met/Trp single-codon families excluded:', sorted(SINGLE), '=>', sorted(CODON_TABLE[c] for c in SINGLE))
print()

# ---- index case, per gene and concatenated
idx = {}
for seg in 'LMS':
    recs = rd(f'{ST}/06_epuyen_orfs/by_patient/Patient_1_{seg}.fasta')
    idx[seg] = recs[0][1]
print('TABLE 2 / TABLE 1 -- ANDV index case (Patient 1)')
per_inc = [cai(idx[s], False) for s in 'LMS']
per_exc = [cai(idx[s], True) for s in 'LMS']
for s, a, b in zip('LMS', per_inc, per_exc):
    print(f'    gene {s}: inc={a:.4f}  exc={b:.4f}  d={b-a:+.4f}')
row('index-case genes (mean of 3)', per_inc, per_exc, 0.6619)
cat = idx['L'] + idx['M'] + idx['S']
row('index case L+M+S (concatenated)', [cai(cat, False)], [cai(cat, True)], 0.647)

# ---- successors
print('\nTABLE 1 -- successor genomes (L+M+S concatenated per patient)')
pats = {}
for f in glob.glob(f'{ST}/06_epuyen_orfs/by_patient/Patient_*_*.fasta'):
    bn = os.path.basename(f)[:-6].split('_')
    pid, seg = int(bn[1]), bn[2]
    pats.setdefault(pid, {})[seg] = rd(f)[0][1]
succ = [p for p in sorted(pats) if p != 1 and set(pats[p]) == {'L','M','S'}]
si = [cai(pats[p]['L']+pats[p]['M']+pats[p]['S'], False) for p in succ]
se = [cai(pats[p]['L']+pats[p]['M']+pats[p]['S'], True) for p in succ]
row(f'successors (n={len(succ)} complete)', si, se, 0.648)
print(f'    inc range {min(si):.4f}-{max(si):.4f}   exc range {min(se):.4f}-{max(se):.4f}')
print(f'    (published range 0.6468-0.6739)')

# ---- panels
print('\nTABLE 1 / TABLE 2 -- human reference panels')
for nm, f, pub in (('lung HKG panel', 'lung_hkg_cds.fasta', 0.7853),
                   ('kidney HKG panel', 'kidney_hkg_cds.fasta', 0.7720)):
    recs = rd(f'{ST}/10_cai_and_coadaptation_results/{f}')
    vi = [cai(s, False) for _, s in recs]
    ve = [cai(s, True) for _, s in recs]
    row(f'{nm} (n={len(recs)})', vi, ve, pub)
    print(f'    inc range {min(vi):.4f}-{max(vi):.4f}   exc range {min(ve):.4f}-{max(ve):.4f}')

# ---- 22 species nucleoprotein (Table 7)
print('\nTABLE 7 -- 22-species nucleoprotein')
recs = rd(f'{ST}/04_genus_alignments/S_segment/genus_diversity_S_cds.fasta')
vi = [cai(s, False) for _, s in recs]
ve = [cai(s, True) for _, s in recs]
row(f'all species (n={len(recs)})', vi, ve)
print(f'    max |delta| over species: {max(abs(a-b) for a,b in zip(vi,ve)):.4f}')
print(f'\n  per-species (header, inc, exc, delta):')
for (h, _), a, b in sorted(zip(recs, vi, ve), key=lambda t: t[0][0]):
    print(f'    {h[:46]:<46} {a:.4f}  {b:.4f}  {b-a:+.4f}')
