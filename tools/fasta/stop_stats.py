#!/usr/bin/env python
import argparse
from Bio import SeqIO

STOP_CODON_NAMES = {"TAG": "Amber", "TAA": "Ochre", "TGA": "Opal"}


def extract_stops(fasta):
    codon_usage = {}
    for record in SeqIO.parse(fasta, "fasta"):
        seq = str(record.seq)[-3:]
        try:
            codon_usage[seq] += 1
        except KeyError:
            codon_usage[seq] = 1

    for key in sorted(codon_usage):
        yield (STOP_CODON_NAMES.get(key.upper(), "None"), key, str(codon_usage[key]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Summarise stop codon usage", epilog=""
    )
    parser.add_argument(
        "fasta", type=argparse.FileType("r"), help="CDS Sequence (DNA + stop)"
    )
    args = parser.parse_args()
    print("# Name\tCodon\tCount")
    for (name, key, value) in extract_stops(**vars(args)):
        print("{}\t{}\t{}".format(name, key, value))
