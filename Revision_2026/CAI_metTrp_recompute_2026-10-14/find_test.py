import glob, os, sys
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'cai_recompute.py')).read().split("print('Met/Trp")[0])
from scipy import stats
from statistics import mean
PATH={'HTNV','SEOV','DOBV','PUUV','SNV','ANDV','BCCV','BAYV','CHOV','LANV'}
recs=rd(f'{ST}/04_genus_alignments/S_segment/genus_diversity_S_cds.fasta')
g={h.split('|')[0]:(cai(s,False),cai(s,True)) for h,s in recs}
for mode,i in (('INC',0),('EXC',1)):
    a=[v[i] for k,v in g.items() if k in PATH]
    b=[v[i] for k,v in g.items() if k not in PATH]
    print(f"--- {mode}  path mean={mean(a):.4f} env mean={mean(b):.4f}")
    cands={
      'MWU two-sided':stats.mannwhitneyu(a,b,alternative='two-sided').pvalue,
      'MWU greater':stats.mannwhitneyu(a,b,alternative='greater').pvalue,
      'MWU 2s no-correct':stats.mannwhitneyu(a,b,alternative='two-sided',use_continuity=False).pvalue,
      'Kruskal-Wallis':stats.kruskal(a,b).pvalue,
      't-test (equal var)':stats.ttest_ind(a,b).pvalue,
      't-test (Welch)':stats.ttest_ind(a,b,equal_var=False).pvalue,
      't-test one-sided':stats.ttest_ind(a,b,alternative='greater').pvalue,
    }
    for k,v in cands.items():
        flag=' <== 0.083' if abs(v-0.083)<0.0015 else ''
        print(f"     {k:<22} p={v:.4f}{flag}")
