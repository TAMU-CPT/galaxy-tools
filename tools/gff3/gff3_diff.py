#!/usr/bin/env python
import sys
import logging
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_qual_value
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def gff3_diff(gff3_1, gff3_2):
    recs1 = {}
    recs2 = {}

    for rec1 in GFF.parse(gff3_1):
        for feature in rec1.features:
            if feature.location.strand == 1:
                recs1[feature.location.start] = feature
            else:
                recs1[feature.location.end] = feature

    for rec2 in GFF.parse(gff3_2):
        for feature in rec2.features:
            if feature.location.strand == 1:
                recs2[feature.location.start] = feature
            else:
                recs2[feature.location.end] = feature


    for i in recs1:
        print i, recs2[i]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reports differences between two gff3 files')
    parser.add_argument('gff3_1', type=argparse.FileType("r"), help='first gff3 file')
    parser.add_argument('gff3_2', type=argparse.FileType("r"), help='first gff3 file')
    args = parser.parse_args()
    gff3_diff(**vars(args))
