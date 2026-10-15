"""Emit the revised Supplementary Tables S1 (per-patient) and S7 (per-gene panel).

Two corrections are applied together, because both land in the same cells:
  (1) CAI with Met and Trp excluded (Sharp & Li 1987);
  (2) the lung panel's B2M record harmonised from XM_005254549.3 (a PREDICTED model
      RefSeq, 369 bp) to NM_004048.4 (curated, 360 bp) -- which is what the project's
      own accession_resolution.tsv resolved B2M to, and what the kidney panel already
      used. The same gene entered the two panels as two different records.

Nothing is emitted until the published values are reproduced in include-mode with the
ORIGINAL accessions. ENC is reproduced too, because one ENC cell (lung B2M) changes
with the accession; Met/Trp cannot move ENC, as single-codon families enter none of
its degeneracy classes.

Run:  python3 -I emit_supplementary.py
"""
import sys, os, glob, json
from math import log, exp
from statistics import mean

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# cai_recompute.py prints its own validation report at import time; silence it
# so this script's output stands alone. Importing it rather than copying it is
# deliberate: the CAI and weight code is the already-validated one, unchanged.
import io as _io, contextlib as _ctx
with _ctx.redirect_stdout(_io.StringIO()):
    from cai_recompute import (W, SINGLE, CODON_TABLE, AA_TO_CODONS,
                               cai, strip_stop, rd, ST)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')

# ---------------------------------------------------------------- ENC (Wright 1990)
# Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6, F per family = (n*sum(p^2) - 1)/(n - 1)
DEGEN = {}
for aa, cs in AA_TO_CODONS.items():
    DEGEN.setdefault(len(cs), []).append(aa)


def enc(seq, undefined='as_published'):
    """Wright (1990) Nc.

    `undefined` selects what to do with a degeneracy class whose F cannot be
    estimated (no family with n > 1) or estimates to exactly 0 (a family used
    perfectly evenly, which makes 1/F a division by zero):

      'as_published' : substitute F = 1. This is the original pipeline's
                       convention and reproduces all 30 published panel ENC
                       values, NDUFA3 (38.842) included. Note what it means:
                       F = 1 is MAXIMAL bias, imposed on a family whose codons
                       are maximally even. It is reproduced here rather than
                       corrected, because ENC is independent of the Met/Trp
                       exclusion, so the only ENC cell this revision changes is
                       the lung B2M row, and that cell must be computed on the
                       same convention as the 29 beside it.
      'wright'       : F3 = (F2 + F4)/2, Wright's own suggestion for an absent
                       or singleton Ile family. Reported alongside for NDUFA3.
    """
    seq = strip_stop(seq)
    counts = {}
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        aa = CODON_TABLE.get(c)
        if aa is None or aa == '*':
            continue
        counts.setdefault(aa, {}).setdefault(c, 0)
        counts[aa][c] += 1

    fbar = {}
    for k in (2, 3, 4, 6):
        fs = []
        for aa in DEGEN.get(k, []):
            fam = counts.get(aa)
            if not fam:
                continue
            n = sum(fam.values())
            if n <= 1:
                continue
            s = sum((c / n) ** 2 for c in fam.values())
            fs.append((n * s - 1.0) / (n - 1.0))
        fbar[k] = mean(fs) if fs else None

    # Note `is None or == 0` rather than a falsy test: F == 0.0 is a real,
    # computed estimate (perfectly even usage), not a missing one, and it is
    # the case that actually arises -- in NDUFA3's Ile family.
    unusable = [k for k in (2, 3, 4, 6) if fbar[k] is None or fbar[k] == 0.0]
    if undefined == 'as_published':
        for k in unusable:
            fbar[k] = 1.0
    elif undefined == 'wright':
        if 3 in unusable and fbar[2] and fbar[4]:
            fbar[3] = (fbar[2] + fbar[4]) / 2.0
        for k in (2, 4, 6):
            if fbar[k] is None or fbar[k] == 0.0:
                fbar[k] = 1.0
    nc = 2.0
    for k, w in ((2, 9), (3, 1), (4, 5), (6, 3)):
        if not fbar[k]:
            return float('nan')
        nc += w / fbar[k]
    return min(nc, 61.0)


def rdd(p):
    return dict(rd(p))


# ---------------------------------------------------------------- validation gate
def load_pub(p):
    rows = {}
    with open(p) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            if not line.strip():
                continue
            d = dict(zip(hdr, line.rstrip('\n').split('\t')))
            rows[d['gene_symbol']] = d
    return rows

PUB = {
    'lung': load_pub(f'{ST}/10_cai_and_coadaptation_results/hkg_metrics.tsv'),
    'kidney': load_pub(f'{ST}/10_cai_and_coadaptation_results/kidney_hkg_metrics.tsv'),
}
SEQ = {
    'lung': rdd(f'{ST}/10_cai_and_coadaptation_results/lung_hkg_cds.fasta'),
    'kidney': rdd(f'{ST}/10_cai_and_coadaptation_results/kidney_hkg_cds.fasta'),
}

print('=' * 78)
print('VALIDATION GATE -- reproduce the published panel CAI and ENC, original accessions')
print('=' * 78)
bad = 0
for panel in ('lung', 'kidney'):
    print(f'\n[{panel}]  {len(PUB[panel])} genes')
    print(f'  {"gene":<10} {"accession":<20} {"CAI pub":>8} {"CAI inc":>8}  '
          f'{"ENC pub":>8} {"ENC calc":>8}')
    for hdr, s in SEQ[panel].items():
        sym, acc = hdr.split('_', 1)
        p = PUB[panel].get(sym)
        if p is None:
            print(f'  {sym:<10} {acc:<20} -- NOT IN PUBLISHED TABLE'); bad += 1; continue
        ci, ep = cai(s, False), enc(s)
        dc = abs(ci - float(p['CAI']))
        de = abs(ep - float(p['ENC']))
        flag = ''
        if dc >= 0.00006:
            flag += f'  <-- CAI off {ci-float(p["CAI"]):+.5f}'; bad += 1
        if de >= 0.0006:
            flag += f'  <-- ENC off {ep-float(p["ENC"]):+.4f}'; bad += 1
        print(f'  {sym:<10} {acc:<20} {p["CAI"]:>8} {ci:>8.4f}  '
              f'{p["ENC"]:>8} {ep:>8.3f}{flag}')

print(f'\nGATE: {"PASS" if bad == 0 else "FAIL (%d mismatches)" % bad}')

# The one gene whose ENC depends on the undefined-class convention
_n = SEQ['lung']['NDUFA3_NM_004542.4']
print(f'\n  ENC convention, the one gene it bites (NDUFA3, {len(strip_stop(_n))//3} codons):')
print(f'    published / as_published (F=1 substitution) : {enc(_n, "as_published"):.3f}')
print(f'    Wright 1990 (F3 = (F2+F4)/2)                : {enc(_n, "wright"):.3f}')
print('    Its only 3-fold family (Ile) has 2 codons used once each, so F3 = 0')
print('    exactly and 1/F3 is undefined. NOT changed here: ENC does not depend')
print('    on the Met/Trp exclusion, so this cell is unaffected by the revision.')
if bad:
    print('\nRefusing to emit. The implementation does not reproduce the published table,\n'
          'so no revised value from it can be trusted.')
    sys.exit(1)

# ---------------------------------------------------------------- Table S7
B2M_NEW = 'B2M_NM_004048.4'
print('\n' + '=' * 78)
print('TABLE S7 (revised) -- per-gene panel CAI with Met/Trp excluded, B2M harmonised')
print('=' * 78)

s7 = []
for panel in ('lung', 'kidney'):
    for hdr, s in SEQ[panel].items():
        sym, acc = hdr.split('_', 1)
        harmon = ''
        if panel == 'lung' and sym == 'B2M':
            s = SEQ['kidney'][B2M_NEW]
            acc_new = B2M_NEW.split('_', 1)[1]
            harmon = f'harmonised from {acc} (predicted model RefSeq) to {acc_new}'
            acc = acc_new
        pub = PUB[panel][sym]
        s7.append(dict(
            panel=panel, gene=sym, accession=acc, length_nt=len(strip_stop(s)),
            cai_published=float(pub['CAI']), cai_revised=round(cai(s, True), 4),
            enc_published=float(pub['ENC']), enc_revised=round(enc(s), 3),
            note=harmon))

hdr7 = ['panel', 'gene', 'accession', 'length_nt', 'cai_published', 'cai_revised',
        'enc_published', 'enc_revised', 'note']
p7 = os.path.join(OUT, 'TableS7_panel_cai_revised.tsv')
with open(p7, 'w') as fh:
    fh.write('\t'.join(hdr7) + '\n')
    for r in s7:
        fh.write('\t'.join(str(r[k]) for k in hdr7) + '\n')

for r in s7:
    d = r['cai_revised'] - r['cai_published']
    de = r['enc_revised'] - r['enc_published']
    mark = '  *** ' + r['note'] if r['note'] else ''
    print(f'  {r["panel"]:<7} {r["gene"]:<10} {r["accession"]:<16} '
          f'CAI {r["cai_published"]:.4f} -> {r["cai_revised"]:.4f} ({d:+.4f})   '
          f'ENC {r["enc_published"]:7.3f} -> {r["enc_revised"]:7.3f} ({de:+.3f}){mark}')

for panel in ('lung', 'kidney'):
    v = [r['cai_revised'] for r in s7 if r['panel'] == panel]
    o = [r['cai_published'] for r in s7 if r['panel'] == panel]
    print(f'\n  {panel:<7} n={len(v)}  CAI mean {mean(o):.4f} -> {mean(v):.4f}   '
          f'range {min(o):.4f}-{max(o):.4f} -> {min(v):.4f}-{max(v):.4f}')

# ---------------------------------------------------------------- Table S1
print('\n' + '=' * 78)
print('TABLE S1 (revised) -- per-patient CAI, Met/Trp excluded')
print('=' * 78)

pats = {}
for f in glob.glob(f'{ST}/06_epuyen_orfs/by_patient/Patient_*_*.fasta'):
    bn = os.path.basename(f)[:-6].split('_')
    pats.setdefault(int(bn[1]), {})[bn[2]] = rd(f)[0][1]

s1 = []
for pid in sorted(pats):
    segs = pats[pid]
    present = ''.join(s for s in 'LMS' if s in segs)
    cat = ''.join(segs[s] for s in 'LMS' if s in segs)
    s1.append(dict(
        patient=pid,
        role='index case' if pid == 1 else 'successor',
        segments=present,
        complete='yes' if present == 'LMS' else 'NO',
        length_nt=len(strip_stop(cat)),
        cai_published=round(cai(cat, False), 4),
        cai_revised=round(cai(cat, True), 4),
        cai_L=round(cai(segs['L'], True), 4) if 'L' in segs else 'NA',
        cai_M=round(cai(segs['M'], True), 4) if 'M' in segs else 'NA',
        cai_S=round(cai(segs['S'], True), 4) if 'S' in segs else 'NA'))

hdr1 = ['patient', 'role', 'segments', 'complete', 'length_nt',
        'cai_published', 'cai_revised', 'cai_L', 'cai_M', 'cai_S']
p1 = os.path.join(OUT, 'TableS1_patient_cai_revised.tsv')
with open(p1, 'w') as fh:
    fh.write('\t'.join(hdr1) + '\n')
    for r in s1:
        fh.write('\t'.join(str(r[k]) for k in hdr1) + '\n')

for r in s1:
    w = '' if r['complete'] == 'yes' else '   <-- INCOMPLETE, segments ' + r['segments']
    print(f'  Patient {r["patient"]:>2}  {r["role"]:<11} {r["segments"]:<3} '
          f'{r["length_nt"]:>5} nt  {r["cai_published"]:.4f} -> {r["cai_revised"]:.4f}{w}')

comp = [r for r in s1 if r['complete'] == 'yes' and r['role'] == 'successor']
inc = [r for r in s1 if r['complete'] != 'yes']
cv = [r['cai_revised'] for r in comp]
cp = [r['cai_published'] for r in comp]
print(f'\n  28 patients: {len(comp)} complete successors, 1 index case, '
      f'{len(inc)} incomplete ({", ".join("P%d (%s)" % (r["patient"], r["segments"]) for r in inc)})')
print(f'  complete successors (n={len(comp)}):')
print(f'      published  mean {mean(cp):.4f}  range {min(cp):.4f}-{max(cp):.4f}')
print(f'      revised    mean {mean(cv):.4f}  range {min(cv):.4f}-{max(cv):.4f}')
for r in inc:
    print(f'  Patient {r["patient"]} ({r["segments"]} only): '
          f'published {r["cai_published"]:.4f}, revised {r["cai_revised"]:.4f}')

allsucc = [r for r in s1 if r['role'] == 'successor']
print(f'\n  as-published pooling (all {len(allsucc)} successors, incomplete included):')
print(f'      published  mean {mean([r["cai_published"] for r in allsucc]):.4f}  '
      f'range {min(r["cai_published"] for r in allsucc):.4f}-'
      f'{max(r["cai_published"] for r in allsucc):.4f}')
print(f'      revised    mean {mean([r["cai_revised"] for r in allsucc]):.4f}  '
      f'range {min(r["cai_revised"] for r in allsucc):.4f}-'
      f'{max(r["cai_revised"] for r in allsucc):.4f}')

print(f'\nWROTE\n  {os.path.abspath(p1)}\n  {os.path.abspath(p7)}')
