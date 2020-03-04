#!/usr/bin/env python
import sys
import copy
import argparse
import logging
from Bio import SeqIO

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def split(fasta_file=None, id="merged"):
    for rec in SeqIO.parse(fasta_file, "fasta"):
        for header in rec.description.split(">"):
            nrec = copy.deepcopy(rec)
            nrec.id = header
            nrec.description = ""
            yield nrec


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_file", type=argparse.FileType("r"))

    args = parser.parse_args()
    for seq in split(**vars(args)):
        SeqIO.write([seq], sys.stdout, "fasta")
