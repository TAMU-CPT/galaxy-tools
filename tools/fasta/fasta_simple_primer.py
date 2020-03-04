#!/usr/bin/env python
import argparse
from Bio import SeqIO

import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def generate_primers(fasta_file, end_overlap=500):
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = record.seq
        output = ">%s %s\n%s[%s]%s\n" % (
            record.id,
            record.description,
            seq[-end_overlap:],
            seq[-3:],
            seq[0:end_overlap],
        )
        yield output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simple primer generator")
    parser.add_argument("fasta_file", type=argparse.FileType("r"), help="Fasta file")
    parser.add_argument("--end_overlap", type=int, help="End overlap", default=500)

    args = parser.parse_args()
    for seq in generate_primers(**vars(args)):
        print(seq)
