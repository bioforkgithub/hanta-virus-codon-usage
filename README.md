# Codon usage, host-tRNA compatibility and purifying selection in Andes virus and the genus *Orthohantavirus*

Analysis code for a study of codon usage bias, human tRNA compatibility, RNA
folding stability and selective constraint in hantaviruses, examined at two very
different timescales: within a single contact-traced Andes virus outbreak, and
across 22 named *Orthohantavirus* species.

This repository contains **code only**. No sequence data, results tables,
figures or manuscript files are included — every input is downloaded from public
databases by the scripts themselves (see [Data](#7-data)).

---

## 1. Disclosure: use of artificial intelligence

**Artificial intelligence was used in producing this work, and this section
states exactly how.**

The analysis pipeline in this repository was written with the assistance of
**Claude (Anthropic)**, used through Claude Code as a programming assistant.
That assistance covered:

- drafting and refactoring the Python and R scripts collected here;
- implementing published algorithms in Python — the codon adaptation index
  (Sharp & Li 1987), the effective number of codons (Wright 1990) and the
  Nei–Gojobori (1986) counting method — after the intended EMBOSS binaries on
  the host system were found to be unreliable (see
  [`shared_utils.py`](#shared_utilspy) for that decision and its reasoning);
- organising intermediate results and generating the figure scripts;
- reviewing, correcting and copy-editing the accompanying manuscript, verifying
  every reference DOI against Crossref, and auditing the code against the
  manuscript's own Methods section — a check that found and corrected several
  discrepancies between what the manuscript claimed and what the code does.

**What AI did not do.** No sequence data was generated, altered or imputed by
AI. No result, statistic, table or figure value in this repository or the
associated manuscript was produced by a language model — every number comes from
running these scripts on public database records. No references were invented;
each one was verified against Crossref before use.

**Accountability.** The author has reviewed the code and retains full
responsibility for its correctness, for the analyses it performs, and for all
conclusions drawn from it. AI assistance does not transfer any part of that
responsibility.

This disclosure is made because reproducibility depends on knowing how research
software was built, and because the analysis relies on in-house implementations
of published methods rather than on established packages. Readers are encouraged
to check those implementations — see
[Validating the implementations](#8-validating-the-implementations).

---

## 2. Summary of the work

Hantaviruses (family *Hantaviridae*, order *Elliovirales*, class
*Bunyaviricetes*) are tri-segmented negative-sense RNA viruses carried by
rodents, shrews, moles and bats. They cause two serious human illnesses:
haemorrhagic fever with renal syndrome (HFRS) and hantavirus
(cardio)pulmonary syndrome (HPS). Like all RNA viruses they depend entirely on
the host's translational machinery, so how well their codon usage matches the
host's preferred codons and tRNA pool may shape how efficiently they replicate.

The study asks one question at two timescales:

> **How strongly does purifying selection constrain the hantavirus nucleoprotein
> gene — and does that constraint track a lineage's ability to cause human
> disease?**

**Short timescale — about eleven weeks.** The 2018–2019 Épuyén outbreak in
Chubut Province, Argentina, is the clearest documented instance of
person-to-person hantavirus transmission. Twenty-eight patients have genome
sequences in GenBank — all three segments for all but one, whose L segment was
never deposited. Comparing every successor genome with the index case gives
within-outbreak polymorphism (pn/ps) — selection caught in the act.

**Long timescale — tens of millions of years.** Twenty-two named
*Orthohantavirus* species, split into ten confirmed human pathogens and twelve
found only in wildlife, give 231 pairwise nucleoprotein dN/dS comparisons —
selection long since completed.

Supporting measures computed throughout: **CAI** (codon adaptation index),
**ENC** (effective number of codons), **co-AI** (correlation between codon usage
and human tRNA gene copy number), **mFE** (windowed minimum RNA folding energy),
and synonymous-codon randomisation. Viral genes are benchmarked against two
independently published human housekeeping-gene panels — one from lung, one from
kidney — because Andes virus affects both organs.

### Provenance

The pipeline is a **methodological replication** of

> Victor MP, Das R, Ghosh TC (2023) *Revisiting SARS-CoV-2: evolution,
> polymorphism and compatibility with human tRNA pool in a lineage over a month.*
> Indian Journal of Applied Microbiology 25:11–36

applied to hantaviruses instead of SARS-CoV-2, and extended with a genus-wide
dN/dS comparison and a second (kidney) tissue benchmark that have no counterpart
in the original.

---

## 3. Repository layout

```
.
├── Epuyen_outbreak/              short-timescale analysis (one outbreak)
│   ├── download_epuyen.py
│   ├── download_lung_hkg.py
│   ├── download_kidney_hkg.py
│   ├── build_trna_table.py
│   ├── extract_orfs.py
│   ├── run_all.sh                runs pipeline steps 01–07 in order
│   └── pipeline/
│       ├── shared_utils.py
│       ├── 01_trna_near_cognate.py
│       ├── 02_codon_usage_metrics.py
│       ├── 03_coai.py
│       ├── 04_folding_energy.py
│       ├── 05_randomization.py
│       ├── 06_pnps_mutations.py
│       ├── 07_stats_figures.R
│       ├── 08_kidney_hkg_metrics.py
│       └── 09_tissue_comparison_stats_figures.R
│
├── Genus_Diversity/              long-timescale analysis (22 species)
│   ├── download_genus_diversity.py
│   └── pipeline/
│       ├── 01_translate_align.py
│       ├── 02_pairwise_dnds.py
│       ├── 03_genus_stats_figures.R
│       └── 04_phylogeny.py
│
└── helpers/                      exploratory scripts from earlier stages
    ├── download_sequences.py
    ├── making_one_fasta.py
    ├── filtering_fasta.py
    ├── FOr_Codon_Usage.py
    ├── making_human_CT.py
    └── cai_fasta.py
```

Scripts create their own output directories (`gbk/`, `orfs/`, `results/`,
`figures/`) relative to their parent folder. Those are not tracked here.

---

## 4. Requirements

Two environments were used. The exact versions the published analysis ran on:

**Python 3.14**

| package | version | used for |
|---|---|---|
| biopython | 1.86 | GenBank parsing, FASTA I/O, neighbour-joining tree |
| seqfold | 0.10.2 | RNA minimum free energy folding |
| scipy | 1.17.1 | Pearson correlation, paired *t*-test |
| numpy | 2.4.3 | array support |
| pandas | — | only `helpers/download_sequences.py` |

**R 4.5.3**

| package | version |
|---|---|
| ggplot2 | 4.0.3 |
| ggpubr | 1.0.0 |
| patchwork | 1.3.2 |
| dplyr | 1.2.1 |
| readr | 2.2.0 |
| tidyr | 1.3.2 |
| stringr | 1.6.0 |

**External binary**

| tool | version | used for |
|---|---|---|
| MAFFT | 7.526 | protein alignment of the 22 species |

Quick setup with conda/mamba:

```bash
mamba create -n hanta python=3.14 biopython scipy numpy pandas mafft \
      r-base r-ggplot2 r-ggpubr r-patchwork r-dplyr r-readr r-tidyr r-stringr
mamba activate hanta
pip install seqfold
```

### Configuration

Nothing needs editing. Three optional environment variables let you point the
scripts at your own tools and contact address:

| variable | used by | default |
|---|---|---|
| `ENTREZ_EMAIL` | the five `download_*.py` scripts | the author's address |
| `MAFFT` | `Genus_Diversity/pipeline/01_translate_align.py` | `mafft` on `PATH` |
| `CODON_PY` | `Epuyen_outbreak/run_all.sh` | `python3` on `PATH` |
| `HANTA_R` | `Epuyen_outbreak/run_all.sh` | `Rscript` on `PATH` |

NCBI asks that requests carry a contact address, so please set your own before
running any download script:

```bash
export ENTREZ_EMAIL="you@example.org"
```

If your MAFFT or R live outside `PATH`, point at them directly:

```bash
export MAFFT=/path/to/mafft
export HANTA_R=/path/to/Rscript
```

---

## 5. Program descriptions

### `Epuyen_outbreak/` — the outbreak analysis

#### Data acquisition

**`download_epuyen.py`** — Downloads the Épuyén 2018–2019 Andes virus outbreak
records from GenBank (accessions MN258162–MN258256): each patient's L (RdRp),
M (glycoprotein precursor) and S (nucleoprotein) segment, where deposited. Writes one GenBank
flat file per accession to `gbk/`, plus `manifest.tsv` mapping patient →
segment → accession → collection date. Patient 1 (sampled 2018-11-03) is the
index case; the rest are successors.

**`download_lung_hkg.py`** — Fetches the twelve human lung housekeeping genes
used as the pulmonary benchmark panel, by exact RefSeq accession, and extracts
each coding sequence.

**`download_kidney_hkg.py`** — Fetches the eighteen-gene human kidney reference
panel of Strauss et al. (2021, PLoS ONE 16:e0259373) by resolving each gene
symbol to its canonical RefSeq mRNA, then extracting the single CDS. Andes virus
causes both pulmonary and renal disease, so both tissue panels are relevant.

**`build_trna_table.py`** — Builds a per-codon *Homo sapiens* cytosolic tRNA
gene copy-number table from the GtRNAdb hg38 high-confidence set. Maps anticodon
to codon by reverse complement. Excludes selenocysteine, three entries GtRNAdb
itself flags as anomalous, and counts initiator Met separately from elongator
Met.

**`extract_orfs.py`** — Extracts clean, in-frame ORF nucleotide sequences from
the downloaded GenBank records using the single annotated CDS per segment.
Writes three views — per segment, per patient, and per patient with all ORFs
concatenated — plus `orf_summary.tsv` with length and stop-codon sanity checks.

#### Pipeline

**<a name="shared_utilspy"></a>`pipeline/shared_utils.py`** — Shared library:
FASTA I/O, the codon table, amino-acid polarity classes, CAI, ENC, windowed mFE
and the synonymous-codon randomiser. CAI and ENC are implemented directly in
Python rather than calling EMBOSS, because on the original host `cai` was a
symlink to codonW (a different, incompatible program) and `chips` was a broken
wrapper. This is the most important file to read if you want to check the
methods.

**`pipeline/01_trna_near_cognate.py`** — Fifteen of the 61 sense codons have no
dedicated human tRNA gene and must be read by wobble pairing. This maps each to
its most abundant compatible near-cognate tRNA, following dos Reis et al. (2004).
The mapping is a property of the human tRNA pool, not of any virus.

**`pipeline/02_codon_usage_metrics.py`** — Computes CAI and ENC for each
patient's combined genome, for each individual ORF, and for the lung panel;
writes the per-codon usage tables that step 03 needs.

**`pipeline/03_coai.py`** — The tRNA co-adaptation index: Pearson correlation
between codon usage frequency and near-cognate-mapped human tRNA gene copy
number across all 61 sense codons.

**`pipeline/04_folding_energy.py`** — Windowed minimum folding energy via
seqfold: 50-nucleotide windows, 30-nucleotide step, 37 °C, averaged across each
transcript.

**`pipeline/05_randomization.py`** — Generates 100 synonymous-codon-randomised
versions of each viral ORF, sampling without replacement within each amino-acid
family so that both the protein sequence and the exact codon usage are
preserved — only codon *position* changes. Recomputes mFE for each and compares
with the native sequence. Any difference therefore comes from nucleotide
arrangement alone, not from codon choice.

**`pipeline/06_pnps_mutations.py`** — Nei–Gojobori pn/ps for every
successor-versus-index-case comparison, with pathway averaging for codons
differing at more than one position. Also classifies each substitution by amino
acid change, side-chain polarity change and the tRNA copy-number consequence.
Segments with no deposited sequence are skipped, so the output row count is the
authoritative number of evaluable comparisons.

**`pipeline/07_stats_figures.R`** — *t*-tests, ANOVA and correlations for the
outbreak-level results; produces the CAI, mFE and pn/ps figures.

**`pipeline/08_kidney_hkg_metrics.py`** — Runs the same CAI, ENC, mFE and co-AI
computations over the kidney panel, so the virus can be benchmarked against both
tissues on identical footing.

**`pipeline/09_tissue_comparison_stats_figures.R`** — The three-way comparison:
Andes virus index-case ORFs versus lung panel versus kidney panel, for CAI and
mFE.

### `Genus_Diversity/` — the genus-wide analysis

**`download_genus_diversity.py`** — Downloads one representative RefSeq
S-segment (nucleoprotein) coding sequence for each of 22 curated
*Orthohantavirus* species, tagged pathogenic or environmental. Accessions were
resolved live against RefSeq and manually reviewed for correct segment and
species assignment.

**`pipeline/01_translate_align.py`** — Translates each coding sequence to
protein, aligns the proteins with MAFFT (`--auto`), then back-translates the
alignment to a codon-level nucleotide alignment: each aligned residue maps to
its original codon and gaps become gap codons. The back-translation is
implemented here directly rather than via an external tool such as PAL2NAL.
Aligning on protein and projecting back is more reliable than aligning
nucleotides directly at this depth of divergence.

**`pipeline/02_pairwise_dnds.py`** — Pairwise dN/dS across all
C(22,2) = 231 species pairs, reusing the site- and difference-counting functions
validated on the outbreak data. Columns gapped in either sequence of a pair are
excluded. Each pair is classified pathogenic–pathogenic,
environmental–environmental, or cross.

**`pipeline/03_genus_stats_figures.R`** — Kruskal–Wallis and post-hoc Wilcoxon
tests (Benjamini–Hochberg corrected) comparing dN/dS across the three pair
categories, plus per-species CAI and ENC by clinical category. **Note the caveat
written into the script:** the 231 pairwise values are not independent — each
species appears in 21 pairs — so these are exploratory statistics, not a
phylogenetically corrected test.

**`pipeline/04_phylogeny.py`** — Neighbour-joining tree of the 22 species from
the protein alignment, using uncorrected p-distance, midpoint-rooted. Presented
to illustrate clustering topology only; it is not a molecular-clock phylogeny
and no reported result depends on its branch lengths.

### `helpers/` — earlier exploratory scripts

Kept for completeness. These predate the two pipelines above and are **not**
part of the published analysis; the pipelines do not import them.

| script | what it does |
|---|---|
| `download_sequences.py` | Bulk GenBank download driven by a TSV accession list |
| `making_one_fasta.py` | Concatenates FASTA files found under a directory tree |
| `filtering_fasta.py` | Selects records from a large FASTA by Ensembl ID list |
| `FOr_Codon_Usage.py` | Filters a FASTA to records whose length is a multiple of three |
| `making_human_CT.py` | Parses a pasted human codon usage table into machine-readable form |
| `cai_fasta.py` | Standalone CAI calculator with hardcoded human weights |

---

## 6. How to run

Clone, then make the configuration edits in
[section 4](#configuration).

### The outbreak analysis

```bash
cd Epuyen_outbreak

# one-time data acquisition
python download_epuyen.py
python download_lung_hkg.py
python download_kidney_hkg.py
python build_trna_table.py          # needs GtRNAdb hg38 set in hg38-tRNAs/
python extract_orfs.py

# steps 01-07 in order
bash run_all.sh

# the kidney panel and the three-way tissue comparison
python  pipeline/08_kidney_hkg_metrics.py
Rscript pipeline/09_tissue_comparison_stats_figures.R
```

`run_all.sh` is idempotent — safe to re-run from scratch. Individual steps can
also be run on their own, in numerical order, as long as earlier outputs exist.

### The genus-wide analysis

```bash
cd Genus_Diversity

python  download_genus_diversity.py
python  pipeline/01_translate_align.py     # requires MAFFT
python  pipeline/02_pairwise_dnds.py
Rscript pipeline/03_genus_stats_figures.R
python  pipeline/04_phylogeny.py
```

There is no runner script here; run the four steps in order.

### Outputs

Tab-separated results land in `results/`, figures in `figures/` (PNG at 150 dpi,
plus PDF and 300 dpi TIFF for publication) within each analysis directory.

### Runtime and network

Everything runs on a laptop in a few minutes. The download scripts make many
NCBI Entrez calls and rate-limit themselves — do not run them in parallel, and
expect the first acquisition pass to be the slowest part.

---

## 7. Data

No data is included. Every input is public and is fetched by the download
scripts:

| input | source |
|---|---|
| Épuyén outbreak genomes | GenBank MN258162–MN258256 |
| Lung housekeeping panel (12 genes) | RefSeq, by accession |
| Kidney housekeeping panel (18 genes) | RefSeq, via Strauss et al. 2021 gene list |
| Human tRNA gene copy numbers | GtRNAdb 2.0, hg38 high-confidence set |
| 22 *Orthohantavirus* nucleoprotein CDS | RefSeq / annotated GenBank |

The GtRNAdb hg38 set must be downloaded manually from
<https://gtrnadb.ucsc.edu/> into `Epuyen_outbreak/hg38-tRNAs/` before running
`build_trna_table.py`; the others are automated.

---

## 8. Validating the implementations

CAI, ENC and Nei–Gojobori dN/dS are implemented in this repository rather than
taken from an established package, for the reason given under
[`shared_utils.py`](#shared_utilspy). If you intend to rely on these numbers,
cross-check a subset against independent tools — EMBOSS `cai` and `chips` (a
genuine installation), codonW, or MEGA for dN/dS. Discrepancies are worth
reporting as issues.

---

## 9. Citation

If you use this code, please cite the manuscript (details to be added on
acceptance) and this repository.

Please also cite the methods it implements:

- Sharp PM, Li WH (1987) *Nucleic Acids Res* **15**:1281–1295 — CAI
- Wright F (1990) *Gene* **87**:23–29 — ENC
- Nei M, Gojobori T (1986) *Mol Biol Evol* **3**:418–426 — dN/dS counting
- dos Reis M, Savva R, Wernisch L (2004) *Nucleic Acids Res* **32**:5036–5044 — near-cognate decoding
- Katoh K, Standley DM (2013) *Mol Biol Evol* **30**:772–780 — MAFFT
- Zuker M, Stiegler P (1981) *Nucleic Acids Res* **9**:133–148 — folding algorithm
- Chan PP, Lowe TM (2016) *Nucleic Acids Res* **44**:D184–D189 — GtRNAdb
- Cock PJA et al. (2009) *Bioinformatics* **25**:1422–1423 — Biopython
- Victor MP, Das R, Ghosh TC (2023) *Indian J Appl Microbiol* **25**:11–36 — the framework replicated here

---

## 10. License

Released under the **MIT License** — see [`LICENSE`](LICENSE).

You may use, modify and redistribute this code freely, including commercially,
provided the copyright notice and licence text are retained. The code is
provided without warranty.

Note that the licence covers this software only. The sequence records it
downloads are subject to the terms of NCBI, RefSeq and GtRNAdb.
