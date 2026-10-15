#!/usr/bin/env python3
"""
Phase 2 step 2: codon-aware (translate - align - back-translate) alignment of the
within-species nucleoprotein CDS sets.

Same strategy as Work1/Genus_Diversity/pipeline/01_translate_align.py: the CDS is
translated, the PROTEINS are aligned with mafft, and the protein alignment is
back-translated onto the original codons, so no indel can ever shift the frame.

Checks enforced per species (fatal unless listed as a warning):
  - every input CDS length is a multiple of 3
  - no internal stop codon (terminal stop is stripped before translation)
  - the back-translated alignment length is a multiple of 3
  - each row's ungapped codon count equals its input codon count

mafft invocation is fixed (FFT-NS-i, --retree 2 --maxiterate 1000) rather than
--auto so that the algorithm does not change with n; mafft is deterministic for
fixed input and these settings.

Usage:  python3 -I 11_align.py <ANALYSIS_ROOT> <MAFFT_BIN>
Deterministic: no RNG.
"""
import sys, os, subprocess, json

STOPS = {"TAA", "TAG", "TGA"}
CODON_TABLE = {}
_bases = "TCAG"
_aas = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
_i = 0
for b1 in _bases:
    for b2 in _bases:
        for b3 in _bases:
            CODON_TABLE[b1 + b2 + b3] = _aas[_i]
            _i += 1


def read_fasta(path):
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def translate(nt):
    aa = []
    for i in range(0, len(nt), 3):
        c = nt[i:i + 3]
        aa.append(CODON_TABLE.get(c, "X"))
    return "".join(aa)


def main():
    root, mafft = sys.argv[1], sys.argv[2]
    setdir = os.path.join(root, "results", "sets")
    aldir = os.path.join(root, "results", "alignments")
    os.makedirs(aldir, exist_ok=True)
    os.chmod(aldir, 0o700)
    report = {}

    for fn in sorted(os.listdir(setdir)):
        if not fn.endswith(".lato.fasta"):
            continue
        abbrev = fn.split(".")[0]
        recs = read_fasta(os.path.join(setdir, fn))
        if len(recs) < 2:
            report[abbrev] = dict(n=len(recs), skipped="fewer than 2 sequences")
            print("%-5s SKIP n=%d" % (abbrev, len(recs)))
            continue
        warn = []
        prots, codons = [], {}
        for name, nt in recs:
            if len(nt) % 3:
                raise SystemExit("%s %s CDS length %d not a multiple of 3"
                                 % (abbrev, name, len(nt)))
            aa = translate(nt)
            if aa.endswith("*"):
                aa, nt_use = aa[:-1], nt[:-3]
            else:
                nt_use = nt
                warn.append("%s has no terminal stop" % name)
            if "*" in aa:
                raise SystemExit("%s %s has an internal stop codon at aa %d"
                                 % (abbrev, name, aa.index("*") + 1))
            prots.append((name, aa))
            codons[name] = [nt_use[i:i + 3] for i in range(0, len(nt_use), 3)]

        pin = os.path.join(aldir, "%s.prot.fasta" % abbrev)
        with open(pin, "w") as fh:
            for name, aa in prots:
                fh.write(">%s\n%s\n" % (name, aa))
        pout = os.path.join(aldir, "%s.prot.aln.fasta" % abbrev)
        cmd = [mafft, "--retree", "2", "--maxiterate", "1000", "--anysymbol",
               "--quiet", "--thread", "16", pin]
        with open(pout, "w") as fh, open(os.path.join(root, "logs",
                                                      "mafft_%s.log" % abbrev), "w") as lg:
            subprocess.run(cmd, stdout=fh, stderr=lg, check=True)

        aln = read_fasta(pout)
        L = len({len(s) for _, s in aln})
        if L != 1:
            raise SystemExit("%s ragged protein alignment" % abbrev)
        nt_out = []
        for name, aa in aln:
            cs, k, row = codons[name], 0, []
            for ch in aa:
                if ch == "-":
                    row.append("---")
                else:
                    row.append(cs[k]); k += 1
            if k != len(cs):
                raise SystemExit("%s %s back-translation used %d of %d codons"
                                 % (abbrev, name, k, len(cs)))
            nt_out.append((name, "".join(row)))
        alen = len(nt_out[0][1])
        if alen % 3:
            raise SystemExit("%s codon alignment length %d not a multiple of 3"
                             % (abbrev, alen))
        with open(os.path.join(aldir, "%s.codon.aln.fasta" % abbrev), "w") as fh:
            for name, s in nt_out:
                fh.write(">%s\n%s\n" % (name, s))
        report[abbrev] = dict(n=len(recs), aln_nt=alen, aln_codons=alen // 3,
                              input_codon_lengths=sorted({len(v) for v in codons.values()}),
                              warnings=warn, mafft_cmd=" ".join(cmd[:-1]) + " <prot>")
        print("%-5s n=%-4d codons=%-5d input codon lengths=%s %s"
              % (abbrev, len(recs), alen // 3,
                 sorted({len(v) for v in codons.values()}),
                 ("WARN %d" % len(warn)) if warn else ""))

    with open(os.path.join(aldir, "align_report.json"), "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
