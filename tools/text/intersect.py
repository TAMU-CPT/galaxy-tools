#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def joinStripRow(row):
    return '\t'.join(row).strip()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('a', type=argparse.FileType("r"))
    parser.add_argument('a_col', type=int, default=1)
    parser.add_argument('b', type=argparse.FileType("r"))
    parser.add_argument('b_col', type=int, default=1)
    args = parser.parse_args()

    data_a = [
        x.split('\t')
        for x in args.a.readlines()
    ]

    data_b = [
        x.split('\t')
        for x in args.b.readlines()
    ]

    data_a_indexed = {
        row[args.a_col - 1]: row
        for row in data_a
    }

    for row in data_b:
        row_id = row[args.b_col - 1]
        if row_id in data_a_indexed:
            print joinStripRow(data_a_indexed[row_id]), joinStripRow(row)
