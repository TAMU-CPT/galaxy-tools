#!/usr/bin/env python
import sys
import logging
import copy
import argparse
from BCBio import GFF
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def fixed_feature(rec):
    for idx, feature in enumerate(rec.features):
        old_quals = feature.qualifiers
        feature.qualifiers = {
            'ID': ['biopy-%s' % idx],
            'product': ['tRNA-' + feature.qualifiers['Codon'][0]],
            'note': ['anticodon: %s' % old_quals['Anticodon'][0].upper()],
            'experiment': old_quals['source'],
        }

        gene = copy.deepcopy(feature)
        gene.qualifiers = {}
        gene.type = 'gene'
        gene.qualifiers['ID'] = 'gene-%s' % idx
        gene.qualifiers['source'] = 'aragorn'

        gene.sub_features = [feature]
        yield gene


def gff_filter(gff3):
    for rec in GFF.parse(gff3):
        rec.features = sorted(list(fixed_feature(rec)), key=lambda x: x.location.start)
        rec.annotations = {}
        GFF.write([rec], sys.stdout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='add parent gene features to CDSs')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='GFF3 annotations')
    args = parser.parse_args()
    gff_filter(**vars(args))
