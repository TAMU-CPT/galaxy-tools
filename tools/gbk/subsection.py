#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def genbank_subsection(genbank_files, table=None, coordinates=None, include_unlisted=False):
    if (table is None and coordinates is None) or \
            (table is not None and coordinates is not None):
        raise Exception("Must specify a coordinate table (X)OR enter coordinates manually")

    if table is not None:
        cut_sites = {}
        for row in table.readlines():
            (a, b, c) = row.split('\t')

            b = int(b)
            c = int(c.strip())
            cut_sites[a] = sorted([b, c])
    else:
        start, end = map(int, coordinates.split(','))

    for genbank_file in genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            if table is not None:
                # If found, cut
                rid = record.id

                if rid in cut_sites:
                    start, end = cut_sites[rid]
                    yield [record[start - 1:end]]
                # if want unlisted, print
                elif include_unlisted:
                    yield [record]
            else:
                yield [record[start - 1:end]]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract subsets of genbank files', epilog='')
    parser.add_argument('genbank_files', nargs='+', type=file, help='Genbank files')
    parser.add_argument('--table', type=file, help='Table of coordinates to cut')
    parser.add_argument('--coordinates', help='Manually entered coordinates')
    parser.add_argument('--include_unlisted', action='store_true',
                        help='If coordinates aren\'t listed in the '
                        'file, still include in output')

    args = parser.parse_args()
    for genbank_record in genbank_subsection(**vars(args)):
        SeqIO.write(genbank_record, sys.stdout, 'genbank')
