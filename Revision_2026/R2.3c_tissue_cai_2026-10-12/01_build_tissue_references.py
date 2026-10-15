#!/usr/bin/env python3
"""R2.3c + R4.1 step 1 - derive tissue-specific CAI weight tables.

Satisfies R4.1 (weights derived from highly expressed genes of the tissue,
not a generic human codon table) and supplies R2.3c with lung, kidney AND a
panel of 49 control tissues.

Inputs
  HPA consensus tissue RNA (Gene, Gene name, Tissue, nTPM), downloaded this run
  Ensembl GRCh38 CDS fasta (Work1, read-only) - gene:ENSG... in the header

Method
  per tissue: rank genes by consensus nTPM, keep the top N whose gene has a
  usable CDS, sum codon counts over those CDS, convert to per-amino-acid
  fractions, then w = f / f_max (Sharp & Li). Writes one table per tissue.
"""
import collections
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cai_lib  # noqa: E402

HPA = sys.argv[1]
CDSFA = sys.argv[2]
OUTDIR = sys.argv[3]
TOPN = int(sys.argv[4]) if len(sys.argv) > 4 else 500

# ---- 1. one CDS per gene: longest protein_coding CDS, must be a clean ORF ---
gene_re = re.compile(r"gene:(ENSG\d+)")
bio_re = re.compile(r"transcript_biotype:(\S+)")
sym_re = re.compile(r"gene_symbol:(\S+)")

best = {}
name = None
buf = []
hdr = None


def consider(hdr, seq):
    g = gene_re.search(hdr)
    b = bio_re.search(hdr)
    if not g or not b or b.group(1) != "protein_coding":
        return
    seq = seq.upper()
    if len(seq) < 300 or len(seq) % 3 or set(seq) - set("ACGT"):
        return
    if not seq.startswith("ATG") or seq[-3:] not in cai_lib.STOPS:
        return
    gid = g.group(1)
    cur = best.get(gid)
    if cur is None or len(seq) > len(cur[1]):
        sym = sym_re.search(hdr)
        best[gid] = (sym.group(1) if sym else "", seq)


with open(CDSFA) as fh:
    for line in fh:
        if line.startswith(">"):
            if hdr is not None:
                consider(hdr, "".join(buf))
            hdr = line
            buf = []
        else:
            buf.append(line.strip())
if hdr is not None:
    consider(hdr, "".join(buf))

sys.stderr.write("genes with usable CDS: %d\n" % len(best))

# precompute codon counts per gene once
gene_counts = {g: cai_lib.codon_counts(s) for g, (_sym, s) in best.items()}

# ---- 2. HPA expression ------------------------------------------------------
expr = collections.defaultdict(dict)  # tissue -> gene -> nTPM
with open(HPA) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    i_g, i_t, i_v = head.index("Gene"), head.index("Tissue"), head.index("nTPM")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(i_g, i_t, i_v):
            continue
        try:
            expr[p[i_t]][p[i_g]] = float(p[i_v])
        except ValueError:
            continue

tissues = sorted(expr)
sys.stderr.write("tissues: %d\n" % len(tissues))

os.makedirs(os.path.join(OUTDIR, "tissue_weights"), exist_ok=True)
manifest = {"topN": TOPN, "n_genes_with_cds": len(best), "tissues": {},
            "hpa_file": HPA, "cds_file": CDSFA}

for t in tissues:
    ranked = sorted((v, g) for g, v in expr[t].items() if g in gene_counts)
    ranked.reverse()
    chosen = ranked[:TOPN]
    tot = collections.Counter()
    for _v, g in chosen:
        tot.update(gene_counts[g])
    frac = cai_lib.fractions_from_counts(tot)
    w = cai_lib.weights_proper(frac)
    slug = re.sub(r"[^a-z0-9]+", "_", t.lower()).strip("_")
    path = os.path.join(OUTDIR, "tissue_weights", "w_%s.tsv" % slug)
    with open(path, "w") as out:
        out.write("codon\taa\tfraction\tweight\n")
        for c in sorted(cai_lib.CODON2AA):
            out.write("%s\t%s\t%.6f\t%s\n" % (
                c, cai_lib.CODON2AA[c], frac.get(c, 0.0),
                ("%.6f" % w[c]) if c in w else "NA"))
    manifest["tissues"][t] = {
        "slug": slug, "weights_file": path,
        "n_genes_used": len(chosen),
        "min_nTPM_in_topN": chosen[-1][0] if chosen else None,
        "max_nTPM": chosen[0][0] if chosen else None,
        "total_codons": sum(tot.values()),
        "top10_genes": [best[g][0] or g for _v, g in chosen[:10]],
    }

with open(os.path.join(OUTDIR, "tissue_reference_manifest.json"), "w") as fh:
    json.dump(manifest, fh, indent=2, sort_keys=True)

print("built %d tissue weight tables (topN=%d)" % (len(tissues), TOPN))
for t in ("lung", "kidney"):
    if t in manifest["tissues"]:
        m = manifest["tissues"][t]
        print("%-8s codons=%d  nTPM range %.1f-%.1f  top genes: %s" % (
            t, m["total_codons"], m["min_nTPM_in_topN"], m["max_nTPM"],
            ", ".join(m["top10_genes"][:6])))
