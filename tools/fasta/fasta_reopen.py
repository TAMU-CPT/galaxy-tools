#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def reopen_contig(fasta_file, after):
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = record.seq
        record.seq = seq[after - 1:] + seq[0:after - 1]
        yield [record]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reopen contig')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), help='Fasta file')
    parser.add_argument('after', type=int, help='Reopen contig after this base (1-indexed)')

    args = parser.parse_args()
    for seq in reopen_contig(**vars(args)):
        SeqIO.write(seq, sys.stdout, 'fasta')
