#!/usr/bin/env python
import sys
import logging
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_contains
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def gff_reopen(gff3, index=1):
    for rec in GFF.parse(gff3):
        # Reopen
        if len(list(feature_lambda(rec.features, feature_test_contains, {'index': index}, subfeatures=False))) > 0:
            log.warn("WARNING: Index chosen is in the middle of a feature. This feature will disappear from the output")
        log.debug(rec.annotations)
        # TODO: This call removes metadata!
        rec = rec[index:] + rec[0:index]
        log.debug(rec.annotations)
        yield rec

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reopen a set of GFF3 annotations')
    parser.add_argument('gff3', type=file, help='GFF3 annotations')
    parser.add_argument('index', type=int, help='Index to reopen genome at')
    args = parser.parse_args()

    for rec in gff_reopen(**vars(args)):
        GFF.write([rec], sys.stdout)
