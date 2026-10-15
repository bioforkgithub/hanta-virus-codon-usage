#!/usr/bin/env python3
"""
Panel entry for group `pathogenic-B`, to be added to PANEL in
01_inventory_N_cds.py.

WHY THIS IS A SIDECAR AND NOT AN EDIT.  Two subagents were given different
species groups and both wrote to scripts/01_inventory_N_cds.py.  The
environmental-A version won and is the one to keep: it identifies a complete N
CDS by a PROTEIN-length band plus /coded_by back-mapping, and it records that
title-keyword matching was MEASURED to miss records.  That is a better method
than the nucleotide-CDS-feature matcher the pathogenic-B version used, and it
has actually been run.  So pathogenic-B should be a new PANEL entry in that
script, not a second script.  This file is kept separate so that adding it
cannot clobber work still in progress.

TO USE: paste the dict entry below into PANEL in 01_inventory_N_cds.py, then

  python3 -I 01_inventory_N_cds.py --outdir .../inputs --results .../results \
      --logs .../logs --group pathogenic-B --email <ncbi contact email>

ONE THING TO CHECK FIRST.  The AA_MIN, AA_MAX = 410, 450 band was justified on
the note that "all six species inventoried here translate to exactly 433 aa".
That is true of the environmental panel, not of this one: Andes virus N is
**428 aa** (REVISED_THESIS.md claim 1 counts S as 428 aa across the Epuyen
chain).  428 is inside 410-450, so the band holds -- but it is no longer a
single-value band, so the log should print the realised aa-length distribution
per species rather than assert 433.

EXPECTED COUNTS, for reconciliation.  A prior harvest for the R2.3b analysis
(analysis/R2.3b_isolate_host_cai_2026-10-12/) already pulled all five of these
species on 2026-10-12.  Its complete-clean-N-CDS counts, which additionally
required a resolvable /host and so are LOWER BOUNDS, were:

    ANDV 125   CHOV 12   BAYV 4   LANV 2   BCCV 0

and the nuccore record ceilings across all three segments were:

    ANDV 909   CHOV 89   LANV 51   BAYV 22   BCCV 8

A protein-database run should return at least the lower bounds.  If it returns
fewer, the protein route is missing records and that must be resolved before any
of these counts are used.  See results/NOTES_inventory_pathogenic-B.md.
"""

PANEL_PATHOGENIC_B = {
    "pathogenic-B": [
        ("ANDV", "Andes orthohantavirus",             "NC_003466.1", "pathogenic", "rodent"),
        ("BCCV", "Black Creek Canal orthohantavirus", "AB689163.1",  "pathogenic", "rodent"),
        ("BAYV", "Bayou orthohantavirus",             "NC_038298.1", "pathogenic", "rodent"),
        ("CHOV", "Choclo orthohantavirus",            "NC_038373.1", "pathogenic", "rodent"),
        ("LANV", "Laguna Negra orthohantavirus",      "NC_038505.1", "pathogenic", "rodent"),
    ],
}

# ANDV CARRIES A TAXONOMIC DECISION the other four do not.  An [Organism] query
# is subtree-inclusive, and only 701 of the 909 Andes records carry the
# "Orthohantavirus andesense" organism name; the other ~208 are sub-lineage
# names (Oran, Bermejo, Lechiguanas, Maciel, Pergamino, Araraquara,
# Juquitiba-like, Maripa, Tunari, Central Plata, Castelo dos Sonhos, Neembucu,
# Jabora).  Pooling them would measure BETWEEN-genotype divergence -- the
# saturated quantity that piN/piS exists to escape.  Recommend restricting ANDV
# to the andesense organism name (sensu stricto, matching labels.tsv
# NC_003466.1) and recording the excluded organism names in the log.
ANDV_SENSU_STRICTO_ORGANISM_NAMES = (
    "orthohantavirus andesense",
    "andes virus",
    "andes orthohantavirus",
)
