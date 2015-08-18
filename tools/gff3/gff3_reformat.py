#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
import logging
logging.basicConfig(level=logging.INFO)

def reformat(data):
    for input_file in data:
        for record in GFF.parse(input_file):
            record.annotations = {}
            GFF.write([record], sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reformat GFF files')
    parser.add_argument('data', type=file, help='Input annotations')
    args = parser.parse_args()
    reformat(**vars(args))
