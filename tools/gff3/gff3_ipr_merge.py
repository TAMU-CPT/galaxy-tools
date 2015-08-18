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


def merge_interpro(gff3, interpro):
    ipr_additions = {}
    badkeys = ('Name', 'ID', 'Target', 'date', 'status', 'signature_desc', 'source', 'md5', 'score')

    for rec in GFF.parse(interpro):
        ipr_additions[rec.id] = {}
        for feature in rec.features:
            quals = feature.qualifiers
            for key in quals:
                if key not in ipr_additions[rec.id]:
                    ipr_additions[rec.id][key] = set()
                for value in quals[key]:
                    ipr_additions[rec.id][key].add(value)

        for key in badkeys:
            if key in ipr_additions[rec.id]:
                del ipr_additions[rec.id][key]

    for rec in GFF.parse(gff3):
        for feature in rec.features:
            if feature.id in ipr_additions:
                for key in ipr_additions[feature.id]:
                    if key not in feature.qualifiers:
                        feature.qualifiers[key] = []

                    feature.qualifiers[key] += list(ipr_additions[feature.id][key])
        GFF.write([rec], sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extract features from a GFF3 file based on ID/qualifiers')
    parser.add_argument('gff3', type=file, help='GFF3 annotations')
    parser.add_argument('interpro', type=file, help='GFF3 annotations')
    args = parser.parse_args()
    merge_interpro(**vars(args))
