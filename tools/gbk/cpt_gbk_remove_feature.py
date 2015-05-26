#!/usr/bin/env python
import argparse
from Bio import SeqIO
import StringIO
import logging
logging.basicConfig(level=logging.INFO)


def remove_feature(genbank_file, feature_type):
    feature_types = feature_type.split(',')
    for record in SeqIO.parse(genbank_file, "genbank"):
        record.features = [x for x in record.features if x.type not in feature_types]

        # Print out the genbank file
        handle = StringIO.StringIO()
        SeqIO.write([record], handle, "genbank")
        print handle.getvalue()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Remove specific features from a Genbank file')
    parser.add_argument('genbank_file', type=file, help='Genbank file')
    parser.add_argument('feature_type', help='Feature type to remove')

    args = parser.parse_args()
    remove_feature(**vars(args))
