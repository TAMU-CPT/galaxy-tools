#!/usr/bin/env python
import argparse
from Bio import SeqIO
from zipfile import ZipFile


def spl(genbank):
    with ZipFile('out.zip', 'w') as zipfile:
        for record in SeqIO.parse(genbank, 'genbank'):
            filename = record.id + '.gbk'
            with open(filename, 'w') as handle:
                SeqIO.write([record], handle, 'genbank')
            zipfile.write(filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='splits multi-genbank file')
    parser.add_argument('genbank', type=str, help='genbank file')

    args = parser.parse_args()
    spl(args.genbank)
