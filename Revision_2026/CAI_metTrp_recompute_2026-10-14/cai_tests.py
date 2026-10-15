import glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'cai_recompute.py')).read().split("print('Met/Trp")[0])
from scipy.stats import mannwhitneyu
from statistics import mean

# --- Patient 2 check (no L segment)
pats={}
for f in glob.glob(f'{ST}/06_epuyen_orfs/by_patient/Patient_*_*.fasta'):
    bn=os.path.basename(f)[:-6].split('_'); pats.setdefault(int(bn[1]),{})[bn[2]]=rd(f)[0][1]
print("Patient 2 segments present:", sorted(pats[2]))
p2=''.join(pats[2][s] for s in 'LMS' if s in pats[2])
print(f"Patient 2 CAI on available segments only: inc={cai(p2,False):.4f}  exc={cai(p2,True):.4f}")
allsucc=[p for p in sorted(pats) if p!=1]
vi=[cai(''.join(pats[p][s] for s in 'LMS' if s in pats[p]),False) for p in allsucc]
ve=[cai(''.join(pats[p][s] for s in 'LMS' if s in pats[p]),True) for p in allsucc]
print(f"ALL {len(allsucc)} successors (Patient 2 on M+S only, as published):")
print(f"   inc mean={mean(vi):.4f} range {min(vi):.4f}-{max(vi):.4f}   (published 0.648, 0.6468-0.6739)")
print(f"   exc mean={mean(ve):.4f} range {min(ve):.4f}-{max(ve):.4f}")

# --- Table 7 pathogenic vs environmental CAI test
PATH={'HTNV','SEOV','DOBV','PUUV','SNV','ANDV','BCCV','BAYV','CHOV','LANV'}
recs=rd(f'{ST}/04_genus_alignments/S_segment/genus_diversity_S_cds.fasta')
g={}
for h,s in recs:
    ab=h.split('|')[0]
    g[ab]=(cai(s,False),cai(s,True))
print(f"\nTABLE 7 CAI: pathogenic (n={len(PATH)}) vs environmental (n={len(g)-len(PATH)})")
for mode,i in (('INCLUDED (published)',0),('EXCLUDED (revised)',1)):
    a=[v[i] for k,v in g.items() if k in PATH]
    b=[v[i] for k,v in g.items() if k not in PATH]
    u,p=mannwhitneyu(a,b,alternative='two-sided')
    print(f"   {mode:<22} path mean={mean(a):.4f}  env mean={mean(b):.4f}  U={u:.1f}  p={p:.4f}")
print("   (manuscript reports CAI p = 0.083)")
