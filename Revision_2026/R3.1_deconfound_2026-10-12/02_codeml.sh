#!/usr/bin/env bash
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
# R3.1 — head-to-head branch models: DISEASE labelling vs HOST labelling on the
# SAME tree and SAME alignment, under both codon-frequency models.
#
# Settings are byte-identical to R2.2/scripts/05_codeml.sh except for the tree,
# so lnL values are directly comparable with the M0 and two_ratio runs there.
set -u
A=/PATH/TO/Project_Hantavirus/Reviews_ans/analysis/R3.1_deconfound_2026-10-12
# PAML 4.10.10 and IQ-TREE 3.1.1 must be on PATH before running this script.

run() {                       # run <name> <treefile-basename> <model> <CodonFreq>
  local name=$1 tree=$2 model=$3 cf=$4
  local wd="$A/results/codeml/$name"
  mkdir -p "$wd" || return 1
  cd "$wd" || return 1
  cp "$A/inputs/S22_codon_aln.phy" seq.phy
  cp "$A/inputs/$tree" tree.nwk
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
  timeout 7200 codeml codeml.ctl > run.log 2>&1
  echo "done $name exit=$? lnL=$(grep -m1 '^lnL' out.txt 2>/dev/null | awk '{print $5}')"
}
export -f run
export A

# ---- R3.1b head-to-head -----------------------------------------------------
# host2  : foreground = the 12 rodent-borne tips       (host labelling)
# three  : 3 classes, pathogenic / rodent-env / rest   (separates the two)
# cladeNW: foreground = New World HPS clade + internals (clade labelling)
# disease2_tulvpath / disease2_tulvown : R3.3c reassignment sensitivity
JOBS=""
for cf in 2 3; do
  sfx=$([ $cf = 2 ] && echo F3x4 || echo F61)
  for t in host2 three cladeNW disease2_tulvpath disease2_tulvown; do
    JOBS="$JOBS ${t}_${sfx}|tree_${t}.nwk|2|$cf"
  done
done

# ---- R3.1c per-species omega (one foreground tip at a time) -----------------
for cf in 2 3; do
  sfx=$([ $cf = 2 ] && echo F3x4 || echo F61)
  for tip in HTNV SEOV DOBV PUUV SNV ANDV BCCV BAYV CHOV LANV PHV TULV \
             TPMV MJNV NVAV ASAV OXBV RKPV SWSV CBNV QZNV JJUV; do
    JOBS="$JOBS pertip_${tip}_${sfx}|tree_pertip_${tip}.nwk|2|$cf"
  done
done

printf '%s\n' $JOBS | tr '|' ' ' | xargs -P 60 -L1 bash -c 'run "$@"' _
echo ALLDONE
