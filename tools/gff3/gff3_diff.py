#!/usr/bin/env python
import logging
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def find_differences(a, b):
    flags = {'strand': False, 'start': False, 'id': False, 'qualifiers': False}
    if a.id != b.id:
        flags['id'] = True

    if a.qualifiers != b.qualifiers:
        flags['qualifiers'] = True

    if a.location.strand == b.location.strand:
        if a.location.strand == 1:
            if a.location.start != b.location.start:
                flags['start'] = True
        else:
            if a.location.end != b.location.end:
                flags['start'] = True
    else:
        flags['strand'] = True

    return flags


def gff3_diff(gff3_1, gff3_2):
    feats1 = {}
    feats2 = {}
    for rec1 in GFF.parse(gff3_1):
        for feat in feature_lambda(rec1.features, feature_test_type, {'type': 'gene'}, subfeatures=True):
            if feat.location.strand == 1:
                feats1[feat.location.start] = feat
            else:
                feats1[feat.location.end] = feat

    for rec2 in GFF.parse(gff3_2):
        for feat in feature_lambda(rec2.features, feature_test_type, {'type': 'gene'}, subfeatures=True):
            if feat.location.strand == 1:
                feats2[feat.location.start] = feat
            else:
                feats2[feat.location.end] = feat

    no_match = []
    flags_list = {}
    for i in feats1:
        try:
            diffs = find_differences(feats1[i], feats2[i])
            # need to somehow check for subfeatures
            del feats2[i]
            for d in diffs:
                if diffs[d]:
                    flags_list[i] = flags  # noqa HXR: Commented out for linting, please remove when ready.
                    break
        except:
            no_match.append(feats1[i])

    print flags_list
    for nm in no_match:
        print nm
    for f in feats2:
        print feats2[f]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reports differences between two gff3 files')
    parser.add_argument('gff3_1', type=argparse.FileType("r"), help='first gff3 file')
    parser.add_argument('gff3_2', type=argparse.FileType("r"), help='first gff3 file')
    args = parser.parse_args()
    gff3_diff(**vars(args))
