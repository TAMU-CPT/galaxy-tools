#!/usr/bin/env python
import sys
import argparse
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def extractCodonsFromDb(target, codondb):
    # Load target data from codondb
    header = None
    lastline = False
    real_target = ":%s:" % target
    for line in codondb:
        if real_target in line:
            header = "" + line
            lastline = True
            continue

        if lastline:
            break

    if header is None:
        print "Not found"
        sys.exit(1)

    spsum_label = (
        "CGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC "
        "UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC "
        "GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU "
        "CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU "
        "AUA AUC AUU AUG UGG UAA UAG UGA"
    )
    spsum = spsum_label.replace("U", "T").split(" ")
    organism_codon_usage = {k: int(v) for (k, v) in zip(spsum, line.strip().split(" "))}
    return organism_codon_usage


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("target", type=str, help="target organism")
    parser.add_argument(
        "--codondb", type=argparse.FileType("r"), help="Average codon database"
    )
    args = parser.parse_args()

    ocu = extractCodonsFromDb(**vars(args))
    print "# Codon\tCount"
    for (key, value) in ocu.iteritems():
        print "%s\t%s" % (key, value)
