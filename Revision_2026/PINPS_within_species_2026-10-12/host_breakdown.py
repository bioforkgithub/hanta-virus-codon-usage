#!/usr/bin/env python3
"""Host / length / reference-presence breakdown of the usable N CDS per species.

Usage: python3 -I host_breakdown.py <INPUTS_DIR> <LABELS_TSV> ABBREV [ABBREV...]
"""
import sys, os, re, collections

HUMAN_RE = re.compile(r"homo sapiens|\bhuman\b|patient", re.I)


def main():
    inputs, labels_tsv = sys.argv[1], sys.argv[2]
    abbrevs = sys.argv[3:]
    refs = {}
    with open(labels_tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > 2:
                refs[p[0]] = p[2]

    print("abbrev\tn_usable\thuman\tanimal\thost_blank\tref_acc\tref_in_usable"
          "\tlen_distribution\ttop_hosts")
    for a in abbrevs:
        tab = os.path.join(inputs, a, "records_table.tsv")
        with open(tab) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(hdr)}
            rows = [l.rstrip("\n").split("\t") for l in fh]
        usable = [r for r in rows if r[idx["usable"]] == "True"]
        human = animal = blank = 0
        hosts = collections.Counter()
        lens = collections.Counter()
        accs = set()
        for r in usable:
            h = r[idx["host"]].strip()
            accs.add(r[idx["accession"]].split(".")[0])
            lens[r[idx["cds_nt"]]] += 1
            if not h:
                blank += 1
            elif HUMAN_RE.search(h):
                human += 1
                hosts[h] += 1
            else:
                animal += 1
                hosts[h] += 1
        ref = refs.get(a, "")
        ref_in = ref.split(".")[0] in accs
        ld = ", ".join("%s:%d" % (k, v) for k, v in sorted(lens.items(),
                       key=lambda kv: -kv[1]))
        th = "; ".join("%s(%d)" % (k, v) for k, v in hosts.most_common(4))
        print("%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s"
              % (a, len(usable), human, animal, blank, ref, ref_in, ld, th))


if __name__ == "__main__":
    main()
