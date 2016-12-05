#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
import sys
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RevCom genbank files')
    parser.add_argument('gbk_files', type=argparse.FileType("r"), nargs='+', help='Genbank files')

    args = parser.parse_args()
    for gbk_file in args.gbk_files:
        for record in SeqIO.parse(gbk_file, 'genbank'):
            SeqIO.write(record.reverse_complement(id=True, name=True, description=True, features=True, annotations=True, letter_annotations=True, dbxrefs=True), sys.stdout, 'genbank')
