#!/usr/bin/env bash
# PUBLICATION NOTE. Absolute paths in this file were rewritten to /PATH/TO/ when the
# repository was published: they pointed into a private working filesystem. Set them
# to wherever you unpacked the data package (see the DOI in the top-level README).
# The code itself is unchanged from the version that produced the reported values.
# Three-category branch model: pathogenic terminals (#1), environmental terminals
# (#2), internal branches (background #0). This is the clean version of the test
# in R4.3: the two-ratio design pooled environmental terminals WITH all
# internal branches, so it confounded pathogenic-vs-environmental with
# terminal-vs-internal (audit finding E1-bis).
set -u
D=/PATH/TO/Project_Hantavirus/Reviews_ans/analysis/R2.2_phylo_dnds_2026-10-12
# PAML 4.10.10 and IQ-TREE 3.1.1 must be on PATH before running this script.
for cf in 2 3; do
  wd="$D/results/codeml_3cat/cf$cf"; mkdir -p "$wd"; cd "$wd" || exit
  cp "$D/inputs/S22_codon_aln.phy" seq.phy; cp "$D/inputs/S22_3cat.nwk" tree.nwk
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
        model = 2
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
  echo "cf$cf exit=$? $(grep -m1 '^lnL' out.txt | sed 's/.*): *\(-*[0-9.]*\).*/lnL=\1/') $(grep -m1 'w (dN/dS) for branches' out.txt)"
done
echo THREECAT_DONE
