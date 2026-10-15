#!/usr/bin/env bash
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
# codeml runs for R2.2 / R3.2 / R4.3: tree-based codon models replacing
# pairwise counting. Each model gets its own directory.
set -u
D=/PATH/TO/Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12
# PAML 4.10.10 and IQ-TREE 3.1.1 must be on PATH before running this script.

run() {                       # run <name> <treefile> <model> <CodonFreq>
  local name=$1 tree=$2 model=$3 cf=$4
  local wd="$D/results/codeml/$name"
  mkdir -p "$wd" && cd "$wd"
  cp "$D/inputs/S22_codon_aln.phy" seq.phy
  cp "$D/inputs/$tree" tree.nwk
  cat > codeml.ctl <<EOF
      seqfile = seq.phy
     treefile = tree.nwk
      outfile = out.txt
        noisy = 3
      verbose = 1
      runmode = 0
      seqtype = 1
    CodonFreq = $cf
        clock = 0
       aaDist = 0
        model = $model
      NSsites = 0
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.4
        getSE = 0
 RateAncestor = 0
   Small_Diff = .5e-6
    cleandata = 0
       method = 0
EOF
  echo "--- $name (model=$model CodonFreq=$cf tree=$tree)"
  timeout 3600 codeml codeml.ctl > run.log 2>&1
  echo "    exit=$? lnL=$(grep -m1 '^lnL' out.txt 2>/dev/null | awk '{print $5}')"
}

# M0: single omega across the whole tree. CodonFreq=2 (F3x4) and 3 (F61)
# -- R3.2 asks specifically for "codeml M0 with a codon-frequency model".
run m0_F3x4     S22_ml_clean.nwk   0 2
run m0_F61      S22_ml_clean.nwk   0 3

# Two-ratio branch model: pathogenic TERMINAL branches as foreground (#1).
# Terminal branches are used because no single bipartition separates the
# pathogenic species (verified), so no clade-level assignment is available.
run two_ratio_F3x4  S22_pathtips.nwk 2 2
run two_ratio_F61   S22_pathtips.nwk 2 3

# Free-ratio: one omega per branch. Descriptive -- heavily parameterised at
# 448 codons, used only to compare terminal-branch omega distributions.
run free_ratio_F3x4 S22_ml_clean.nwk 1 2

echo "done"
