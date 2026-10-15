#!/usr/bin/env bash
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
# R3.1 — convergence check. The F61 fits produced a pattern that cannot come
# from a well-behaved likelihood surface: in the three-class model PHV alone
# gives omega 0.00080 and TULV alone 0.01083, yet the two POOLED give 0.04478,
# which is above both. Re-fit the key F61 models from several starting omegas;
# if the lnL and the omegas move, the F61 results are local optima and must not
# be reported as findings.
set -u
A=/PATH/TO/Project_Hantavirus/Reviews_ans/analysis/R3.1_deconfound_2026-10-12
# PAML 4.10.10 and IQ-TREE 3.1.1 must be on PATH before running this script.
run() {
  local name=$1 tree=$2 model=$3 cf=$4 w0=$5
  local wd="$A/results/codeml_conv/$name"; mkdir -p "$wd"; cd "$wd" || return 1
  cp "$A/inputs/S22_codon_aln.phy" seq.phy; cp "$A/inputs/$tree" tree.nwk
  printf '      seqfile = seq.phy\n     treefile = tree.nwk\n      outfile = out.txt\n        noisy = 3\n      verbose = 1\n      runmode = 0\n      seqtype = 1\n    CodonFreq = %s\n        clock = 0\n       aaDist = 0\n        model = %s\n      NSsites = 0\n        icode = 0\n    fix_kappa = 0\n        kappa = 2\n    fix_omega = 0\n        omega = %s\n        getSE = 0\n RateAncestor = 0\n   Small_Diff = .5e-6\n    cleandata = 0\n       method = 0\n' "$cf" "$model" "$w0" > codeml.ctl
  timeout 10800 codeml codeml.ctl > run.log 2>&1
  echo "done $name lnL=$(grep -m1 '^lnL' out.txt|awk '{print $5}') w=$(grep -m1 'for branches' out.txt|cut -d: -f2)"
}
export -f run; export A
printf '%s\n' \
 "m0_F61_w0.05 S22_ml_clean.nwk 0 3 0.05" \
 "m0_F61_w1.0 S22_ml_clean.nwk 0 3 1.0" \
 "disease2_F61_w0.05 tree_disease2.nwk 2 3 0.05" \
 "disease2_F61_w1.0 tree_disease2.nwk 2 3 1.0" \
 "host2_F61_w0.05 tree_host2.nwk 2 3 0.05" \
 "host2_F61_w1.0 tree_host2.nwk 2 3 1.0" \
 "three_F61_w0.05 tree_three.nwk 2 3 0.05" \
 "three_F61_w1.0 tree_three.nwk 2 3 1.0" \
 "three_phvonly_F61_w0.05 tree_three_phvonly.nwk 2 3 0.05" \
 "cladeNW_F61_w0.05 tree_cladeNW.nwk 2 3 0.05" \
 "disease2_F3x4_w1.0 tree_disease2.nwk 2 2 1.0" \
 "three_F3x4_w1.0 tree_three.nwk 2 2 1.0" \
 | xargs -P 12 -L1 bash -c 'run "$@"' _
echo CONVDONE
