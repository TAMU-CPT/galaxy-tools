#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('a', type=file)
    parser.add_argument('b', type=file)
    args = vars(parser.parse_args())

    data_a = set([
        x.strip()
        for x in args.a.readlines()
    ])

    data_b = set([
        x.strip()
        for x in args.b.readlines()
    ])

    for hit in data_a.intersection(data_b):
        print hit
