#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type, feature_test_true

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    parser.add_argument("changeTo", type=str, help="Feature type to change to")
    parser.add_argument("changeList", type=str, nargs="+", help="Feature types to change")
    args = parser.parse_args()

    
    for record in GFF.parse(args.gff3):
        for feature in feature_lambda(record.features, feature_test_true, {}):
            if feature.type in args.changeList:
                feature.type = args.changeTo
        GFF.write([record], sys.stdout)
