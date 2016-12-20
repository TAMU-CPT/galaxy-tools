#!/usr/bin/env python
import sys
import logging
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def gff3_diff(gff3_1, gff3_2):
    recs1 = {}
    recs2 = {}
    for rec1 in GFF.parse(gff3_1):
        for feat in feature_lambda(rec1.features, feature_test_type, {'type': 'gene'}, subfeatures=True):
            if feat.location.strand == 1:
                recs1[feat.location.start] = feat
            else:
                recs1[feat.location.end] = feat

    for rec2 in GFF.parse(gff3_2):
        for feat in feature_lambda(rec2.features, feature_test_type, {'type': 'gene'}, subfeatures=True):
            if feat.location.strand == 1:
                recs2[feat.location.start] = feat
            else:
                recs2[feat.location.end] = feat


    print len(recs1)
    print len(recs2)
    for i in recs1:
        a = recs1[i]
        try:
            b = recs2[i]
        except:
            print recs1[i]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reports differences between two gff3 files')
    parser.add_argument('gff3_1', type=argparse.FileType("r"), help='first gff3 file')
    parser.add_argument('gff3_2', type=argparse.FileType("r"), help='first gff3 file')
    args = parser.parse_args()
    gff3_diff(**vars(args))
