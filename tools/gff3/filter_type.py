#!/usr/bin/env python
import sys
import argparse
from cpt_gffParser import gffParse, gffWrite
from gff3 import feature_lambda, feature_test_type

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    parser.add_argument("types", type=str, nargs="+", help="Feature type to filter on")
    parser.add_argument("--invert", action="store_true")
    args = parser.parse_args()

    for rec in gffParse(args.gff3):
        tempFeats = feature_lambda(
            rec.features,
            feature_test_type,
            {"types": args.types},
            invert=args.invert,
            subfeatures=False,
        )
        rec.features = []
        for x in tempFeats:
          rec.features.append(x)
        for x in rec.features:
          if "Parent" in x.qualifiers.keys():
            found = 0
            for seek in x.qualifiers["Parent"]:
              for y in rec.features:
                if y.id == seek:
                  found += 1
                  break
            if found < len(x.qualifiers["Parent"]):
              del x.qualifiers["Parent"]
        gffWrite([rec], sys.stdout)
