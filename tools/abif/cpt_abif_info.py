#!/usr/bin/env python
import argparse
from Bio import SeqIO
import logging
logging.basicConfig(level=logging.INFO)


def plot(abif_file):
    records = list(SeqIO.parse(abif_file, "abi"))
    for record in records:
        for key in sorted(record.annotations['abif_raw'].keys()):
            print "%s\t%s" % (key,  record.annotations['abif_raw'][key])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Dump information about AB1 file', epilog='')
    parser.add_argument('abif_file', type=argparse.FileType('rb'),
                        help='ABIF/AB1 file')
    args = parser.parse_args()
    plot(**vars(args))
