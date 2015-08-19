#!/usr/bin/env python
import sys
import argparse
import copy
import logging
logging.basicConfig(level=logging.INFO)
from BCBio import GFF
from Bio import SeqIO
"""
Merge GFF features into a GenBank File
======================================

Useful to add features from GFF producing analysis tools to a GenBank File.

Currently only supports GenBank files with single records (i.e. you cannot
export a GenBank DB, produce gff3 against many subfeatures, and then expect
those to be merged correctly)

"""


def extract_features(gff3_file):
    if gff3_file is None:
        raise ValueError("Must specify gff file")

    for rec in GFF.parse(gff3_file):
        for feat in rec.features:
            if feat.type == 'remark':
                continue

            if feat.type not in ('CDS', 'RBS', "gene", 'terminator'):
                feat.type = 'CDS'
            feat.qualifiers['color'] = ['255 0 0']

            # Remove keys with '-'
            quals = copy.deepcopy(feat.qualifiers)
            for key in quals.keys():
                if '-' in key:
                    del quals[key]
            feat.qualifiers = quals

            yield feat


def merge_features(features=None, genbank_file=None):
    records = SeqIO.parse(genbank_file, 'genbank')

    for record in records:
        for feature in features:
            record.features.append(feature)
        yield [record]


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Merge GFF3 data into a Genbank file')
    parser.add_argument('genbank_file', type=file, help='Genbank file')
    parser.add_argument('gff3_file', type=file, help='GFF3 Input')

    args = parser.parse_args()
    features = extract_features(args.gff3_file)
    for record in merge_features(features=features, genbank_file=args.genbank_file):
        SeqIO.write(record, sys.stdout, "genbank")
