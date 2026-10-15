#!/usr/bin/env bash
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
# R3.1 — does the deconfounded (within-rodent) disease contrast rest on a
# single species? PHV-only and TULV-only versions of the 3-class model.
# TULV-only is already covered by disease2_tulvown / disease2_tulvpath.
set -u
A=/PATH/TO/Project_Hantavirus/Reviews_ans/analysis/R3.1_deconfound_2026-10-12
# PAML 4.10.10 and IQ-TREE 3.1.1 must be on PATH before running this script.
source /dev/stdin <<'FN'
run() {
  local name=$1 tree=$2 model=$3 cf=$4
  local wd="$A/results/codeml/$name"; mkdir -p "$wd"; cd "$wd" || return 1
  cp "$A/inputs/S22_codon_aln.phy" seq.phy; cp "$A/inputs/$tree" tree.nwk
  printf '      seqfile = seq.phy\n     treefile = tree.nwk\n      outfile = out.txt\n        noisy = 3\n      verbose = 1\n      runmode = 0\n      seqtype = 1\n    CodonFreq = %s\n        clock = 0\n       aaDist = 0\n        model = %s\n      NSsites = 0\n        icode = 0\n    fix_kappa = 0\n        kappa = 2\n    fix_omega = 0\n        omega = 0.4\n        getSE = 0\n RateAncestor = 0\n   Small_Diff = .5e-6\n    cleandata = 0\n       method = 0\n' "$cf" "$model" > codeml.ctl
  timeout 7200 codeml codeml.ctl > run.log 2>&1
  echo "done $name lnL=$(grep -m1 '^lnL' out.txt 2>/dev/null | awk '{print $5}')"
}
FN
export -f run; export A
printf '%s\n' \
  "three_phvonly_F3x4 tree_three_phvonly.nwk 2 2" \
  "three_phvonly_F61 tree_three_phvonly.nwk 2 3" \
  "disease2_phvpath_F3x4 tree_disease2_phvpath.nwk 2 2" \
  "disease2_phvpath_F61 tree_disease2_phvpath.nwk 2 3" \
  | xargs -P 4 -L1 bash -c 'run "$@"' _
echo PHVDONE
