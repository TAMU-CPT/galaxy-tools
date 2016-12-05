#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import sys
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def gen_annotation(record, start, end):
    return SeqFeature(
        FeatureLocation(start, end),
        type="source",
        strand=0,
        qualifiers={
            'note': [
                'Source genome originally %s bp' % (len(record.seq)),
                'Cut from %s to %s' % (start + 1, end),
            ]
        }
    )


def genbank_subsection(genbank_files, table=None, coordinates=None, include_unlisted=False):
    if (table is None and coordinates is None) or \
            (table is not None and coordinates is not None):
        raise Exception("Must specify a coordinate table (X)OR enter coordinates manually")

    if table is not None:
        cut_sites = {}
        for row in table.readlines():
            tmp = row.strip().split('\t')

            b = int(tmp[1])
            c = int(tmp[2])

            if tmp[0] not in cut_sites:
                cut_sites[tmp[0]] = []

            cut_sites[tmp[0]].append(sorted([b, c]))
    else:
        start, end = map(int, coordinates.split(','))

    for genbank_file in genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            if table is not None:
                # If found, cut
                rid = record.name

                if rid in cut_sites:
                    for (start, end) in cut_sites[rid]:
                        record.features = [gen_annotation(record, start - 1, end)] + record.features
                        yield [record[start - 1:end]]
                # if want unlisted, print
                elif include_unlisted:
                    record.features = [gen_annotation(record, start - 1, end)] + record.features
                    yield [record]
            else:
                record.features = [gen_annotation(record, start - 1, end)] + record.features
                yield [record[start - 1:end]]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract subsets of genbank files', epilog='')
    parser.add_argument('genbank_files', nargs='+', type=argparse.FileType("r"), help='Genbank files')
    parser.add_argument('--table', type=argparse.FileType("r"), help='Table of coordinates to cut')
    parser.add_argument('--coordinates', help='Manually entered coordinates')
    parser.add_argument('--include_unlisted', action='store_true',
                        help='If coordinates aren\'t listed in the '
                        'file, still include in output')

    args = parser.parse_args()
    for genbank_record in genbank_subsection(**vars(args)):
        SeqIO.write(genbank_record, sys.stdout, 'genbank')
