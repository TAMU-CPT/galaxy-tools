#!/usr/bin/env python
import sys
import logging
import argparse
from BCBio import GFF
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def gff_filter(gff3):
    for rec in GFF.parse(gff3):
        rec.annotations = {}
        GFF.write([rec.reverse_complement(id=True, name=True, description=True, features=True, annotations=True, letter_annotations=True, dbxrefs=True)], sys.stdout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='reverse and complement as set of GFF3 annotations')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='GFF3 annotations')
    args = parser.parse_args()

    gff_filter(**vars(args))
