#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gff3', type=file, help='GFF3 annotations')
    parser.add_argument('types', type=str, nargs='+', help='Feature type to filter on')
    parser.add_argument('--invert', action='store_true')
    args = parser.parse_args()

    for rec in GFF.parse(args.gff3):
        rec.features = feature_lambda(
            rec.features,
            feature_test_type,
            {'types': args.types},
            invert=args.invert,
            subfeatures=False,
        )
        GFF.write([rec], sys.stdout)
