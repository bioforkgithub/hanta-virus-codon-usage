# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'cai_recompute.py')).read().split("print('Met/Trp")[0])
PATH={'HTNV','SEOV','DOBV','PUUV','SNV','ANDV','BCCV','BAYV','CHOV','LANV'}
recs=rd(f'{ST}/04_genus_alignments/S_segment/genus_diversity_S_cds.fasta')
out=['abbrev\tspecies\taccession\tcategory\tCAI_published_incl\tCAI_revised_excl\tdelta']
rows=[]
for h,s in recs:
    p=h.split('|'); ab=p[0]
    i,e=cai(s,False),cai(s,True)
    rows.append((ab,p[1],p[2] if len(p)>2 else '',('pathogenic' if ab in PATH else 'environmental'),i,e))
for r in sorted(rows,key=lambda t:(t[3]!='pathogenic',t[0])):
    out.append(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.4f}\t{r[5]:.4f}\t{r[5]-r[4]:+.4f}")
open('/PATH/TO/Project_Hantavirus/Reviews_ans/analysis/CAI_metTrp_recompute_2026-10-14/results/table7_cai_revised.tsv','w').write('\n'.join(out)+'\n')
print('\n'.join(out))
