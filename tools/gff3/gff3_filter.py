#!/usr/bin/env python
import sys
import logging
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_qual_value
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

__author__ = "Eric Rasche"
__maintainer__ = "Eric Rasche"
__email__ = "esr@tamu.edu"


def gff_filter(gff3, filter_list, attribute_field='ID', subfeatures=True):
    filter_strings = [line.strip() for line in filter_list]
    for rec in GFF.parse(gff3):
        rec.features = feature_lambda(
            rec.features,
            feature_test_qual_value,
            {'qualifier': attribute_field, 'attribute_list': filter_strings},
            subfeatures=subfeatures
        )
        rec.annotations = {}
        GFF.write([rec], sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extract features from a GFF3 file based on ID/qualifiers')
    parser.add_argument('gff3', type=file, help='GFF3 annotations')
    parser.add_argument('filter_list', type=file, help='Child GFF3 annotations to rebase against parent')
    parser.add_argument('--attribute_field', type=str, help='Column 9 Field to search against', default='ID')
    parser.add_argument('--subfeatures', action='store_true', help='Retain subfeature tree of matched features')
    args = parser.parse_args()
    gff_filter(**vars(args))
