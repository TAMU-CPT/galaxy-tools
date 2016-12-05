#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
import sys
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def remove_feature(genbank_files, feature_types):

    for genbank_file in genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            record.features = [x for x in record.features if x.type not in feature_types]

            yield [record]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Remove specific features from a Genbank file')
    parser.add_argument('genbank_files', nargs='+', type=argparse.FileType("r"), help='Genbank files')
    parser.add_argument('--feature_types', nargs='+', help='Feature types to remove')

    args = parser.parse_args()
    for record in remove_feature(**vars(args)):
        SeqIO.write(record, sys.stdout, 'genbank')
