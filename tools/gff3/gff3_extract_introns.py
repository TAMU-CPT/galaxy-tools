#!/usr/bin/env python
import sys
import logging
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def find_introns(gff3):
    for rec in GFF.parse(gff3):
        genes = list(feature_lambda(rec.features, feature_test_type, {'type': 'gene'}, subfeatures=True))
        for gene in genes:
            cdss = list(feature_lambda(gene.sub_features, feature_test_type, {'type': 'CDS'}, subfeatures=False))
            if len(cdss) > 1:
                for cds in sorted(cdss, key=lambda x: x.location.start):
                    print cds.location.start


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extract introns from gene features that have more than one CDS')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='GFF3 annotations')
    args = parser.parse_args()
    find_introns(**vars(args))
