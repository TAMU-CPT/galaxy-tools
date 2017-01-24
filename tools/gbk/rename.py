#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
import sys
from BCBio import GFF
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Renumber genbank files')
    parser.add_argument('data', type=argparse.FileType("r"), help='Annotations')
    parser.add_argument('--filetype', type=str, help='Filetype')
    parser.add_argument('--new_name', type=str, help='New Name')

    args = parser.parse_args()

    # Iterate over our input data
    if args.filetype == 'gff3':
        it = GFF.parse(args.data)
    else:
        it = SeqIO.parse(args.data, args.filetype)

    for idx, record in enumerate(it):
        # Can only handle a single name
        if idx > 1:
            raise Exception("Too many sequences")

        # Update name
        record.id = args.new_name
        record.name = args.new_name

        # Output according to type
        if args.filetype == 'gff3':
            GFF.write([record], sys.stdout)
        else:
            SeqIO.write([record], sys.stdout, args.filetype)
