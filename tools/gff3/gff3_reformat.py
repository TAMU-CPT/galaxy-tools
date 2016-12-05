#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
import logging
logging.basicConfig(level=logging.INFO)


def reformat(data):
    for record in GFF.parse(data):
        record.annotations = {}
        GFF.write([record], sys.stdout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reformat GFF files')
    parser.add_argument('data', type=argparse.FileType("r"), help='Input annotations')
    args = parser.parse_args()
    reformat(**vars(args))
