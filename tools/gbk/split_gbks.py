#!/usr/bin/env python
import argparse
from Bio import SeqIO


def spl(genbank):
    for record in SeqIO.parse(genbank, 'genbank'):
        with open('gbks/' + record.id + '.gbk', 'w') as handle:
            SeqIO.write([record], handle, 'genbank')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='splits multi-genbank file')
    parser.add_argument('genbank', type=str, help='genbank file')

    args = parser.parse_args()
    spl(args.genbank)
