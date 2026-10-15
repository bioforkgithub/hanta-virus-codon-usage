import sys, numpy as np, importlib.util, itertools
spec=importlib.util.spec_from_file_location("m", sys.argv[1]+"/12_pinps.py")
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def pathway_diffs_nostopexcl():
    SD=np.zeros((m.NC+1,m.NC+1)); ND=np.zeros((m.NC+1,m.NC+1))
    for c1 in m.SENSE:
        for c2 in m.SENSE:
            diff=[p for p in range(3) if c1[p]!=c2[p]]
            if not diff: continue
            allp=[]
            for order in itertools.permutations(diff):
                cur,s,n=c1,0.0,0.0
                for p in order:
                    nxt=cur[:p]+c2[p]+cur[p+1:]
                    if m.CODON_AA[nxt]==m.CODON_AA[cur]: s+=1
                    else: n+=1
                    cur=nxt
                allp.append((s,n))
            SD[m.CIDX[c1],m.CIDX[c2]]=sum(p[0] for p in allp)/len(allp)
            ND[m.CIDX[c1],m.CIDX[c2]]=sum(p[1] for p in allp)/len(allp)
    return SD,ND
SS=m.site_counts(False); NS=np.where(np.arange(m.NC+1)<m.NC,3.0-SS,0.0)
SDa,NDa=pathway_diffs_nostopexcl()
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
rng=np.random.default_rng(7); codons=m.SENSE
maxrel=0
for trial in range(25):
    L=200
    a=[codons[rng.integers(0,61)] for _ in range(L)]; b=list(a)
    for _ in range(int(L*0.6)):
        j=int(rng.integers(0,L)); p=int(rng.integers(0,3)); nb="ACGT"[int(rng.integers(0,4))]
        cand=b[j][:p]+nb+b[j][p+1:]
        if m.CODON_AA[cand]!="*": b[j]=cand
    s1="".join(a); s2="".join(b)
    A,NTM=m.encode([s1,s2])
    pN,pS,_,_=m.pairwise(A,NTM,SS,NS,SDa,NDa)
    dN,dS=cal_dn_ds(CodonSeq(s1),CodonSeq(s2),method="NG86")
    jcS=-0.75*np.log(1-4/3*pS[0,1]); jcN=-0.75*np.log(1-4/3*pN[0,1])
    rel=max(abs(jcS-dS)/dS, abs(jcN-dN)/dN)
    maxrel=max(maxrel,rel)
print("25 random 200-codon pairs, NG86 sites=/3 + all-shortest-pathways:")
print("  max relative deviation from Biopython 1.83 cal_dn_ds NG86 = %.3e"%maxrel)
