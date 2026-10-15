#!/usr/bin/env bash
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
# Is the SIGN of the pathogenic-vs-background omega contrast robust to CodonFreq?
# CodonFreq: 0=equal, 1=F1x4, 2=F3x4, 3=F61.  Runs M0 and two-ratio for each.
set -u
D=/PATH/TO/Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12
# PAML 4.10.10 and IQ-TREE 3.1.1 must be on PATH before running this script.
run() {
  local name=$1 tree=$2 model=$3 cf=$4
  local wd="$D/results/codeml_sweep/$name"; mkdir -p "$wd"; cd "$wd" || return
  cp "$D/inputs/S22_codon_aln.phy" seq.phy; cp "$D/inputs/$tree" tree.nwk
  cat > codeml.ctl <<CTL
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
CTL
  timeout 5400 codeml codeml.ctl > run.log 2>&1
  echo "$name done exit=$?"
}
for cf in 0 1; do
  run "m0_cf$cf"        S22_ml_clean.nwk 0 $cf
  run "two_ratio_cf$cf" S22_pathtips.nwk 2 $cf
done
# free-ratio under F61: per-branch omega, less dependent on one fg/bg split
run free_ratio_F61 S22_ml_clean.nwk 1 3
echo SWEEP_COMPLETE
