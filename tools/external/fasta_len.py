#!/usr/bin/env python
import argparse
from Bio import SeqIO
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=argparse.FileType("r"))
    args = parser.parse_args()

    for record in SeqIO.parse(args.fasta, "fasta"):
        sys.stdout.write("%s\t%s\n" % (record.id, len(record)))
