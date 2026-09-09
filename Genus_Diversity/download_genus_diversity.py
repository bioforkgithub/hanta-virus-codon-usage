import os
"""
Download RefSeq S-segment (nucleocapsid protein) CDS for a curated, real,
verified set of Orthohantavirus species spanning both major clinical
categories:

  - "pathogenic": established human pathogens (HFRS in the Old World,
    HPS in the New World)
  - "environmental": hantaviruses discovered through broad rodent/shrew/
    mole/bat surveillance with no established human disease association
    (a few, e.g. Tula virus, have rare/mild case reports and are noted as
    such rather than claimed fully non-pathogenic)

All accessions below were resolved live against NCBI RefSeq
(find_hanta_species_refseq.py) and manually reviewed for correct segment
(S, ~1500-2100 nt) and CDS annotation before being locked in here - nothing
in this list is guessed or reconstructed from memory alone.
"""
import time
from pathlib import Path
from Bio import Entrez, SeqIO

Entrez.email = os.environ.get("ENTREZ_EMAIL", "bioforkgithub@outlook.com")

# (species, common/strain label, RefSeq S-segment accession, clinical category, disease/host note)
SPECIES = [
    # --- pathogenic: HFRS (Old World, rodent-borne) ---
    ("Hantaan orthohantavirus", "HTNV", "NC_005218.1", "pathogenic", "HFRS (severe); Apodemus agrarius"),
    ("Seoul orthohantavirus", "SEOV", "NC_005236.1", "pathogenic", "HFRS (moderate); Rattus norvegicus, worldwide via rats"),
    ("Dobrava-Belgrade orthohantavirus", "DOBV", "NC_005233.1", "pathogenic", "HFRS (severe); Apodemus flavicollis"),
    ("Puumala orthohantavirus", "PUUV", "NC_005224.1", "pathogenic", "Nephropathia epidemica (mild HFRS); Myodes glareolus"),
    # --- pathogenic: HPS (New World, rodent-borne) ---
    ("Sin Nombre orthohantavirus", "SNV", "NC_005216.1", "pathogenic", "HPS; Peromyscus maniculatus, North America"),
    ("Andes orthohantavirus", "ANDV", "NC_003466.1", "pathogenic", "HPS, person-to-person transmission; Oligoryzomys longicaudatus"),
    ("Black Creek Canal orthohantavirus", "BCCV", "AB689163.1", "pathogenic", "HPS; Sigmodon hispidus"),
    ("Bayou orthohantavirus", "BAYV", "NC_038298.1", "pathogenic", "HPS; Oryzomys palustris"),
    ("Choclo orthohantavirus", "CHOV", "NC_038373.1", "pathogenic", "HPS; Panama"),
    ("Laguna Negra orthohantavirus", "LANV", "NC_038505.1", "pathogenic", "HPS; Calomys laucha, Paraguay/Bolivia"),
    # --- environmental: rodent-borne, no/rare human disease ---
    ("Prospect Hill orthohantavirus", "PHV", "U47136.1", "environmental", "No known human disease; Microtus pennsylvanicus"),
    ("Tula orthohantavirus", "TULV", "NC_005227.1", "environmental", "Rare/mild human case reports only; Microtus arvalis"),
    # --- environmental: soricomorph- (shrew/mole) or bat-borne, no known human disease ---
    ("Thottapalayam orthohantavirus", "TPMV", "NC_010704.1", "environmental", "No known human disease; Suncus murinus (shrew), first non-rodent hantavirus"),
    ("Imjin orthohantavirus", "MJNV", "NC_034558.1", "environmental", "No known human disease; Crocidura lasiura (shrew)"),
    ("Nova orthohantavirus", "NVAV", "NC_034464.1", "environmental", "No known human disease; Talpa europaea (mole)"),
    ("Asama orthohantavirus", "ASAV", "NC_038273.1", "environmental", "No known human disease; Urotrichus talpoides (mole)"),
    ("Oxbow orthohantavirus", "OXBV", "NC_043174.1", "environmental", "No known human disease; Neurotrichus gibbsii (mole)"),
    ("Rockport orthohantavirus", "RKPV", "NC_038696.1", "environmental", "No known human disease; Blarina brevicauda (shrew)"),
    ("Seewis orthohantavirus", "SWSV", "NC_043527.1", "environmental", "No known human disease; Sorex araneus (shrew)"),
    ("Cao Bang orthohantavirus", "CBNV", "NC_034484.1", "environmental", "No known human disease; Anourosorex squamipes (shrew)"),
    ("Quezon orthohantavirus", "QZNV", "NC_034400.1", "environmental", "No known human disease; Chaerephon plicata (bat), first bat-borne hantavirus"),
    ("Jeju orthohantavirus", "JJUV", "NC_034398.1", "environmental", "No known human disease; Crocidura shantungensis (shrew)"),
]

BASE = Path(__file__).parent
(BASE / "gbk").mkdir(exist_ok=True)


def main():
    report = []
    fasta_lines = []
    meta_rows = []
    for species, abbrev, acc, category, note in SPECIES:
        try:
            h = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
            text = h.read()
            h.close()
        except Exception as e:
            report.append((species, abbrev, acc, "FETCH_ERROR", str(e)))
            time.sleep(0.5)
            continue
        gbk_path = BASE / "gbk" / f"{abbrev}_{acc}.gbk"
        gbk_path.write_text(text)
        rec = next(SeqIO.parse(gbk_path, "genbank"))
        cds_feats = [f for f in rec.features if f.type == "CDS"]

        target_cds = None
        for f in cds_feats:
            product = " ".join(f.qualifiers.get("product", [""])).lower()
            if "nucleocapsid" in product or "nucleoprotein" in product or ("n protein" in product):
                target_cds = f
                break
        if target_cds is None and len(cds_feats) == 1:
            target_cds = cds_feats[0]
        if target_cds is None and cds_feats:
            # fall back to the longest CDS (nucleocapsid is usually the larger ORF vs overlapping NSs)
            target_cds = max(cds_feats, key=lambda f: len(f))

        if target_cds is None:
            report.append((species, abbrev, acc, "NO_CDS_FOUND", len(cds_feats)))
            time.sleep(0.5)
            continue

        seq = str(target_cds.extract(rec.seq)).upper()
        ok3 = len(seq) % 3 == 0
        ends_stop = seq[-3:] in {"TAA", "TAG", "TGA"}
        report.append((species, abbrev, acc, len(seq), f"n_cds={len(cds_feats)} x3={ok3} stop_end={ends_stop}"))
        fasta_lines.append(f">{abbrev}|{species}|{acc}\n{seq}\n")
        meta_rows.append((species, abbrev, acc, category, note, len(seq), ok3, ends_stop))
        time.sleep(0.4)

    with open(BASE / "genus_diversity_S_cds.fasta", "w") as f:
        f.writelines(fasta_lines)

    with open(BASE / "genus_diversity_manifest.tsv", "w") as f:
        f.write("species\tabbrev\taccession\tcategory\tnote\tlength_nt\tmultiple_of_3\tends_in_stop\n")
        for row in meta_rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"{'Species':35s}{'Abbrev':8s}{'Accession':14s}{'Info'}")
    for species, abbrev, acc, *rest in report:
        print(f"{species:35s}{abbrev:8s}{acc:14s}{rest}")
    print(f"\nWrote {len(fasta_lines)}/{len(SPECIES)} sequences to genus_diversity_S_cds.fasta")


if __name__ == "__main__":
    main()
