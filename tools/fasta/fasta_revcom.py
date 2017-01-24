#!/usr/bin/env python
import argparse
import sys
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RevCom fasta files')
    parser.add_argument('fa_files', type=argparse.FileType("r"), nargs='+', help='fasta files')

    kwargs = dict(
        id=True, name=True, description=True, features=True, annotations=True,
        letter_annotations=True, dbxrefs=True
    )

    args = parser.parse_args()
    for fa_file in args.fa_files:
        for record in SeqIO.parse(fa_file, 'fasta'):
            SeqIO.write(record.reverse_complement(**kwargs), sys.stdout, 'fasta')
