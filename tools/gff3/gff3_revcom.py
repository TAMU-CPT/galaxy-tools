#!/usr/bin/env python
import sys
import logging
import argparse
from BCBio import GFF
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

__author__ = "Eric Rasche"
__maintainer__ = "Eric Rasche"
__email__ = "esr@tamu.edu"


def gff_filter(gff3):
    for rec in GFF.parse(gff3):
        rec.annotations = {}
        GFF.write([rec.reverse_complement(id=rec.id)], sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='reverse and complement as set of GFF3 annotations')
    parser.add_argument('gff3', type=file, help='GFF3 annotations')
    args = parser.parse_args()

    gff_filter(**vars(args))
