#!/usr/bin/env python
import sys
import logging
import argparse
from BCBio import GFF
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

__author__ = "Eric Rasche"
__maintainer__ = "Eric Rasche"
__email__ = "esr@tamu.edu"


def gff_filter(gff3, attribute_field='ID', filter_list):
    filter_strings = [line.strip() for line in filter_list]
    for rec in GFF.parse(parent):
        rec.features = [
            feat for feat in rec.features if
            # if the attribute_value is in the filter_strings list
            any([attribute_value in filter_strings for attribute_value in feat.get(attribute_field, [])])
        ]
        GFF.write([rec], sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extract features from a GFF3 file based on ID/qualifiers')
    parser.add_argument('gff3', type=file, help='GFF3 annotations')
    parser.add_argument('attribute_field', type=str, help='Column 9 Field to search against', default='ID')
    parser.add_argument('filter_list', type=file, help='Child GFF3 annotations to rebase against parent')
    args = parser.parse_args()
    gff_filter(**vars(args))
