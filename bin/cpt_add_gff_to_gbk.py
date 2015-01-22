#!/usr/bin/env python
import argparse
import StringIO
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


def extract_features(gff3_file=None, fasta_file=None, **kwd):
    if gff3_file is None:
        raise ValueError("Must specify gff file")

    features = []
    for rec in GFF.parse(gff3_file):
        features.extend(rec.features)
    return features

def merge_features(features=None, genbank_file=None, **kwd):
    records = SeqIO.parse(genbank_file, 'genbank')
    output = StringIO.StringIO()

    for record in records:
        record.features.extend(features)
        SeqIO.write(record, output, "genbank")

    return output.getvalue()


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Merge GFF3 data into a Genbank file')
    parser.add_argument('genbank_file', type=file, help='Genbank file')
    parser.add_argument('gff3_file', type=file, help='GFF3 Input')
    parser.add_argument('fasta_file', type=file, nargs='?', help='Fasta file used in GFF3 producing analysis')

    args = vars(parser.parse_args())
    features = extract_features(**args)
    result = merge_features(features=features, **args)
    print result
