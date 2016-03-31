#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gff3', type=file, help='GFF3 annotations')

    for rec in GFF.parse(args.gff3):
        rec.annotations = {}

        GFF.write([rec], sys.stdout)
