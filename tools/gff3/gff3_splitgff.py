#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from BCBio import GFF

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('data', type=argparse.FileType("r"), help='GFF3 File')
    parser.add_argument('--gff', type=argparse.FileType('w'), help='Output Annotations', default='data.gff3')
    parser.add_argument('--fasta', type=argparse.FileType('w'), help='Output Sequence', default='data.fa')
    args = parser.parse_args()

    for record in GFF.parse(args.data):
        GFF.write([record], args.gff)
        record.description = ""
        SeqIO.write([record], args.fasta, 'fasta')
        sys.exit()
