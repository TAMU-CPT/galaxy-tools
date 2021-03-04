#!/usr/bin/env python
import sys
import argparse
from cpt_gffParser import gffParse, gffWrite
import logging

logging.basicConfig(level=logging.INFO)


def reformat(data):
    for record in gffParse(data):
        record.annotations = {}
        gffWrite([record], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat GFF files")
    parser.add_argument("data", type=argparse.FileType("r"), help="Input annotations")
    args = parser.parse_args()
    reformat(**vars(args))
