#!/usr/bin/env python
import argparse
import StringIO
from Bio import SeqIO

def genbank_subsection(genbank_file=None, table=None, coordinates=None, include_unlisted=False):
    if (table is None and coordinates is None) or \
            (table is not None and coordinates is not None):
        raise Exception("Must specify a coordinate table (X)OR enter coordinates manually")

    records = SeqIO.parse(genbank_file, "genbank")


    if table is not None:
        cut_sites = {}
        for row in table.readlines():
            (a, b, c) = row.split('\t')
            cut_sites[a] = (int(b), int(c.strip()))

        for record in records:
            output = StringIO.StringIO()
            # If found, cut
            if record.id in cut_sites:
                start, end = cut_sites[record.id]
                SeqIO.write(record[start-1:end], output, "genbank")
            # if want unlisted, print
            elif include_unlisted:
                SeqIO.write(record, output, "genbank")
            #Otherwise, ignore
            print output.getvalue()
    else:
        start, end = map(int, coordinates.split(','))
        for record in records:
            output = StringIO.StringIO()
            SeqIO.write(record[start-1:end], output, "genbank")
            print output.getvalue()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract subsets of genbank files', epilog='')
    parser.add_argument('genbank_file', type=file, help='Genbank file')
    parser.add_argument('--table', help='Table of coordinates to cut')
    parser.add_argument('--coordinates', help='Manually entered coordinates')
    parser.add_argument('--include_unlisted', action='store_true',
                        help='If coordinates aren\'t listed in the file, still include in output')

    args = parser.parse_args()
    shinefind(**vars(args))
