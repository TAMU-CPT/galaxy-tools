#!/usr/bin/env python
from BCBio import GFF
import argparse
import sys


def parse_gff(locations, gff3):
    locs = []
    for line in locations:
        # Consume lines from tabular list of base locations and convert into int list
        line = line.strip()
        if line:
            locs.append(int(line))
    #sort for speed
    locs.sort()

    for rec in GFF.parse(gff3):
        matched_features = []
        for feat in rec.features:
            for loc in locs:
                if loc in feat:
                    # base location is found within this feature's boudary
                    matched_features.append(feat)
                elif loc > feat.location.end:
                    # locations are now beyond this feature, skip checking
                    break
        rec.features = matched_features
        GFF.write([rec], sys.stdout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract features which contain one or more base locations')
    parser.add_argument('locations', type=argparse.FileType("r"),
                        help='table of newline separated base locations')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='GFF3 annotations')
    args = parser.parse_args()
    parse_gff(**vars(args))
