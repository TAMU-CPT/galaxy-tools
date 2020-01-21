#!/usr/bin/env python
import sys
import argparse
import logging
from Bio import SeqIO

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def drop_id(fasta_file=None):
    for rec in SeqIO.parse(fasta_file, "fasta"):
        rec.description = ''
        ind = str(rec.seq).find("##")
        if ind != -1: # This method causes mid-file comments (such as from Apollo sequences) to be appended to the end of the previous sequence
          rec.seq = rec.seq[0:ind]
        yield rec


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), help='Genbank file')

    args = parser.parse_args()
    for rec in drop_id(**vars(args)):
        SeqIO.write([rec], sys.stdout, "fasta")
