import os
"""
Download a Human Kidney Housekeeping Gene (Kidney HKG / K-HKG) panel, by resolving
each gene symbol to its canonical human RefSeq mRNA accession via NCBI, then
extracting each gene's single CDS.

Gene panel source: Strauss P, Mikkelsen H, Furriol J. "Variable expression of
eighteen common housekeeping genes in human non-cancerous kidney biopsies."
PLoS One. 2021;16(12):e0259373. (18 candidate reference genes screened for
expression stability across 30 kidney-tissue expression datasets; YWHAZ and
SLC4A1AP were the two found most stable.)

Andes virus (and hantaviruses generally) cause both pulmonary (HPS) and renal
(HFRS / AKI) disease, so this kidney-tissue-matched panel complements the
Lung HKG panel already downloaded (download_lung_hkg.py).
"""
import time
from pathlib import Path
from Bio import Entrez, SeqIO

Entrez.email = os.environ.get("ENTREZ_EMAIL", "bioforkgithub@outlook.com")

GENE_SYMBOLS = [
    "GAPDH", "ACTG1", "REL", "ACTB", "SLC4A1AP", "YWHAZ", "RPS13", "NOP10",
    "PGAM1", "PPIA", "GUSB", "TBP", "RPL13", "HNRNPL", "PCBP1", "RER1",
    "PLA2G4A", "B2M",
]

BASE = Path(__file__).parent
OUT_DIR = BASE.parent / "Kidney_HKG"
OUT_DIR.mkdir(exist_ok=True)
(OUT_DIR / "gbk").mkdir(exist_ok=True)


def resolve_accession(symbol):
    """Find a canonical human RefSeq mRNA (NM_) accession for a gene symbol."""
    term = f"{symbol}[Gene Name] AND Homo sapiens[Organism] AND srcdb_refseq[Properties] AND biomol_mrna[Properties]"
    h = Entrez.esearch(db="nuccore", term=term, retmax=50, sort="relevance")
    rec = Entrez.read(h)
    h.close()
    ids = rec["IdList"]
    if not ids:
        return None
    h = Entrez.esummary(db="nuccore", id=",".join(ids))
    summaries = Entrez.read(h)
    h.close()
    # Prefer NM_ (curated mRNA), shortest/lowest accession number = often the
    # principal/oldest-established transcript variant.
    nm_candidates = [s for s in summaries if str(s.get("Caption", "")).startswith("NM_")]
    pool = nm_candidates if nm_candidates else summaries
    pool_sorted = sorted(pool, key=lambda s: (str(s.get("Caption", "")), -int(s.get("Length", 0))))
    return str(pool_sorted[0]["Caption"]), str(pool_sorted[0].get("Title", ""))


def main():
    cds_fasta_lines = []
    report = []
    resolved_map = []
    for symbol in GENE_SYMBOLS:
        resolved = resolve_accession(symbol)
        time.sleep(0.4)
        if resolved is None:
            report.append((symbol, None, "NO_HIT", None, None))
            continue
        acc, title = resolved
        resolved_map.append((symbol, acc, title))

        h = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
        text = h.read()
        h.close()
        gbk_path = OUT_DIR / "gbk" / f"{symbol}_{acc}.gbk"
        gbk_path.write_text(text)

        rec = next(SeqIO.parse(gbk_path, "genbank"))
        cds_feats = [f for f in rec.features if f.type == "CDS"]
        if len(cds_feats) != 1:
            report.append((symbol, acc, rec.id, "UNEXPECTED_CDS_COUNT", len(cds_feats)))
            time.sleep(0.4)
            continue
        seq = str(cds_feats[0].extract(rec.seq)).upper()
        ok3 = len(seq) % 3 == 0
        ends_stop = seq[-3:] in {"TAA", "TAG", "TGA"}
        report.append((symbol, acc, rec.id, len(seq), f"x3={ok3} stop_end={ends_stop}"))
        cds_fasta_lines.append(f">{symbol}_{rec.id}\n{seq}\n")
        time.sleep(0.4)

    with open(OUT_DIR / "kidney_hkg_cds.fasta", "w") as f:
        f.writelines(cds_fasta_lines)

    with open(OUT_DIR / "accession_resolution.tsv", "w") as f:
        f.write("symbol\taccession\ttitle\n")
        for symbol, acc, title in resolved_map:
            f.write(f"{symbol}\t{acc}\t{title}\n")

    print(f"{'Symbol':10s}{'Accession':15s}{'Resolved':15s}{'Info'}")
    for symbol, acc, resolved, *rest in report:
        print(f"{symbol:10s}{str(acc):15s}{str(resolved):15s}{rest}")


if __name__ == "__main__":
    main()
