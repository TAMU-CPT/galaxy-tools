#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
import logging
logging.basicConfig(level=logging.INFO)

def reformat(data):
    GFF.write(list(GFF.parse(data)), sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reformat GFF files')
    parser.add_argument('data', type=file, help='Input annotations')
    args = parser.parse_args()
    reformat(**vars(args))
