#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def extract_by_acc(parent=None, acc=None, strict=False):
    acc_list = []
    for line in acc.readlines():
        acc_list.append(line.strip())

    for record in SeqIO.parse(parent, "genbank"):
        acc = record.id

        if strict:
            found = any([x == acc for x in acc_list])
        else:
            found = any([x in acc for x in acc_list])

        if found:
            yield [record]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Subset a genbank file')
    parser.add_argument('parent', type=file, help='Multi-record Genbank file')
    parser.add_argument('acc', type=file, help='Accession list')
    parser.add_argument('--strict', action='store_true', help='Accession numbers '
                        'have versions, setting strict indicates that an exact '
                        'match is required')

    args = parser.parse_args()

    for record in extract_by_acc(**vars(args)):
        SeqIO.write(record, sys.stdout, 'genbank')
