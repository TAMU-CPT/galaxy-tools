#!/usr/bin/env python
import argparse
import logging
from BCBio import GFF
from gff3 import feature_lambda, feature_test_true

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def extract(qualifier, gff3):
    for record in GFF.parse(gff3):
        for feature in feature_lambda(record.features, feature_test_true, {}):
            if qualifier in feature.qualifiers:
                print("%s\t%s" % (feature.id, feature.qualifiers[qualifier][0]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract specified qualifers from features in GFF3", epilog=""
    )
    parser.add_argument("gff3", help="GFF3 File")
    parser.add_argument("qualifier", help="Sepcific qualifier to extract")
    args = parser.parse_args()
    extract(**vars(args))
