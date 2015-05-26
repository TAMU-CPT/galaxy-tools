#!/usr/bin/env python
import argparse
from Bio import SeqIO
import json
import logging
logging.basicConfig(level=logging.INFO)


def plot(abif_file):
    records = list(SeqIO.parse(abif_file, "abi"))
    for record in records:
        print(json.dumps(record.annotations['abif_raw']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Dump information about AB1 file')
    parser.add_argument('abif_file', type=argparse.FileType('rb'),
                        help='ABIF/AB1 file')
    args = parser.parse_args()
    plot(**vars(args))
