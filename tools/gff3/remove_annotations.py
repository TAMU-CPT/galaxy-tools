#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gff3', type=argparse.FileType("r"), help='GFF3 annotations')
    parser.add_argument('--remark', action='store_true', help='Remove remark features')
    args = parser.parse_args()

    for rec in GFF.parse(args.gff3):
        rec.annotations = {}
        if args.remark:
            rec.features = [x for x in rec.features if x.type != 'remark']
        GFF.write([rec], sys.stdout)
