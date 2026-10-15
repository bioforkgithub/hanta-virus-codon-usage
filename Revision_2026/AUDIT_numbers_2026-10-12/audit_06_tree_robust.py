#!/usr/bin/env python3
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
"""C19: would a codon-position-partitioned or codon-substitution-model tree
change the topology materially, and would it change the two claims this analysis
rests on?"""
import re, sys
from pathlib import Path

A = Path("/PATH/TO/"
         "Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12")
R = Path("/PATH/TO/"
         "Project_Hantavirus/Reviews_ans/analysis/AUDIT_numbers_2026-10-12/tree_robust")
PATH10 = {"HTNV", "SEOV", "DOBV", "PUUV", "SNV", "ANDV", "BCCV", "BAYV", "CHOV", "LANV"}
ROD12 = PATH10 | {"PHV", "TULV"}


def parse(nwk):
    s = re.sub(r"\s+", "", nwk.strip().rstrip(";"))
    p = [0]
    def node():
        if s[p[0]] == "(":
            p[0] += 1
            kids = [node()]
            while s[p[0]] == ",":
                p[0] += 1; kids.append(node())
            p[0] += 1
            lbl = ""
            while p[0] < len(s) and s[p[0]] not in "(),":
                lbl += s[p[0]]; p[0] += 1
            return {"kids": kids, "lbl": lbl}
        lbl = ""
        while p[0] < len(s) and s[p[0]] not in "(),":
            lbl += s[p[0]]; p[0] += 1
        return {"kids": [], "lbl": lbl}
    return node()


def leaves(n):
    if not n["kids"]:
        return {n["lbl"].split(":")[0]}
    o = set()
    for k in n["kids"]:
        o |= leaves(k)
    return o


def bips(T, with_support=False):
    allt = leaves(T); out = {}
    def rc(v):
        s = leaves(v)
        if 1 < len(s) < len(allt):
            a, b = frozenset(s), frozenset(allt - s)
            key = a if min(a) < min(b) else b
            m = re.match(r"([\d.]+)/([\d.]+)", v["lbl"])
            out[key] = (float(m.group(1)), float(m.group(2))) if m else None
        for k in v["kids"]:
            rc(k)
    rc(T)
    return out if with_support else set(out)


def clade_possible(T, X):
    allv = leaves(T); ok = [False]
    for Z in (X, allv - X):
        def rc(v):
            lv = leaves(v)
            if Z <= lv:
                if Z == lv:
                    ok[0] = True
                elif v["kids"]:
                    ks = [leaves(k) for k in v["kids"]]
                    if all((s <= Z) or not (s & Z) for s in ks):
                        ok[0] = True
                for k in v["kids"]:
                    rc(k)
        rc(T)
    return ok[0]


ref = parse((A / "results" / "S22_ml.treefile").read_text())
bref = bips(ref)
print("reference: GTR+F+R4 on concatenated nt (the tree used for PAML)\n")

for name, f, desc in [
        ("part3", R / "part3.treefile",
         "codon-position partitioned, ModelFinder per partition, 1000 UFBoot"),
        ("codonmod", R / "codonmod.treefile",
         "codon substitution model GY+F3X4+G4 (-st CODON), 1000 UFBoot")]:
    if not f.exists():
        print(f"{name}: NOT PRESENT (run did not finish) -> [[RESULT: {desc}]]")
        continue
    T = parse(f.read_text())
    b = bips(T)
    rf = len(bref ^ b)
    print(f"--- {name}: {desc}")
    print(f"    splits {len(b)}; Robinson-Foulds symmetric difference vs reference = {rf}"
          f"  ({rf//2} split(s) differ)" if rf else
          f"    splits {len(b)}; IDENTICAL topology to the reference (RF = 0)")
    for s in sorted(bref - b, key=len):
        print(f"      only in reference : {{{','.join(sorted(s))}}}")
    for s in sorted(b - bref, key=len):
        print(f"      only in {name:9s}: {{{','.join(sorted(s))}}}")
    print(f"    pathogenic-10 monophyletic?  {clade_possible(T, PATH10)}   "
          f"(reference: {clade_possible(ref, PATH10)})")
    print(f"    rodent-12   monophyletic?  {clade_possible(T, ROD12)}   "
          f"(reference: {clade_possible(ref, ROD12)})")
    # largest pathogenic-only clade
    best = None
    for s in b:
        for side in (s, frozenset(leaves(T) - s)):
            if side <= PATH10 and (best is None or len(side) > len(best)):
                best = side
    print(f"    largest pathogenic-only clade: {sorted(best)} ({len(best)}/10)")
    print()

# model fit numbers from the partitioned run, for the record
p = R / "part3.iqtree"
if p.exists():
    t = p.read_text()
    for pat in (r"Log-likelihood of the tree: (\S+)", r"Best-fit model.*?: (\S+)"):
        m = re.search(pat, t)
        if m:
            print(f"part3 {pat.split(':')[0][:30]}: {m.group(1)}")
    m = re.findall(r"^\s+\d+\s+(\S+)\s+\d+\s+\d+\s+(\S+)", t, re.M)
    print("part3 per-partition best models:", m[:3])
