#!/usr/bin/env python
import sys
import argparse
import gffutils

"""
Output differences from old filter_type (one that uses bcbio-gff):
    - 'remark' feature is no longer generated
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gff3', type=str, help='GFF3 annotations')
    parser.add_argument('types', type=str, nargs='+', help='Feature type to filter on')
    parser.add_argument('--invert', action='store_true')
    args = parser.parse_args()

    # not sure about merge strategy
    db = gffutils.create_db(args.gff3, 'test.db', merge_strategy="merge", force=True)
    db = gffutils.FeatureDB('test.db')

    for d in db.directives:
        sys.stdout.write('##{0}\n'.format(d))
    for feature in db.all_features():
        if (feature.featuretype in args.types) ^ args.invert:
            sys.stdout.write(str(feature) + '\n')
