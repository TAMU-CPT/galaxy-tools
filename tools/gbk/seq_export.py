#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export sequences from Genbank files')
    parser.add_argument('genbank_files', nargs='+', type=file, help='Genbank file')
    args = parser.parse_args()

    for gbk in args.genbank_files:
        for seq in SeqIO.parse(gbk, 'genbank'):
            SeqIO.write(seq, sys.stdout, 'fasta')
