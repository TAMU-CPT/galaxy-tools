#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def compare_feature_lists(list_a=[], list_b=[]):
    (both, a_only, b_only) = match_feature_lists(list_a=list_a, list_b=list_b)

    data = {
        'PresentInBoth': {
            'header': ['Feature', 'Strand', 'Identical Locations', 'A Start',
                       'A End', 'B Start', 'B End'],
            'data': [],
        },
        'Unique': {
            'header': ['Parent', 'Feature', 'Strand', 'Start', 'End'],
            'data': [],
        }
    }

    for f_a, f_b in both:
        loc_a = '%s..%s' % (f_a.location.start, f_a.location.end)
        loc_b = '%s..%s' % (f_b.location.start, f_b.location.end)
        data['PresentInBoth']['data'].append([
            f_a.id,
            f_a.strand,
            loc_a == loc_b,
            f_a.location.start, f_a.location.end,
            f_b.location.start, f_b.location.end,
        ])

    data['PresentInBoth']['data'].sort(key=lambda x: x[3])

    for f in a_only:
        data['Unique']['data'].append([
            'File 1',
            f.id,
            f.strand,
            f.location.start,
            f.location.end,
        ])

    for f in b_only:
        data['Unique']['data'].append([
            'File 2',
            f.id,
            f.strand,
            f.location.start,
            f.location.end,
        ])

    data['Unique']['data'].sort(key=lambda x: x[3])

    return data


def match_feature_lists(list_a=[], list_b=[]):
    # Try and match up features (end + strand b/c lazy)
    a_only = []
    b_only = []
    both = []

    def mapify_features(feature_list):
        fmap = {}
        for f in feature_list:
            if f.strand == 1:
                fid = '%s.%s' % (f.location.end, f.strand)
            else:
                fid = '%s.%s' % (f.location.start, f.strand)
            # They better not have two features ending in the same location -_-
            fmap[fid] = f
        return fmap

    f_a_map = mapify_features(list_a)
    f_b_map = mapify_features(list_b)

    # Find those keys in common
    both_list = []
    for key in f_a_map:
        if key in f_b_map:
            both_list.append(key)
            both.append([f_a_map[key], f_b_map[key]])
    # Remove those in both from a and b
    for key in both_list:
        del(f_a_map[key])
        del(f_b_map[key])
    # Those left, copy back to original arrays
    for key in f_a_map:
        a_only.append(f_a_map[key])
    for key in f_b_map:
        b_only.append(f_b_map[key])

    return (both, a_only, b_only)


def get_features_from_gbk(gbk_file=None, feature_type='CDS'):
    records = list(SeqIO.parse(gbk_file, "genbank"))
    return [x for x in records[0].features if x.type == feature_type]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare two genbank files')
    parser.add_argument('gbk1', type=argparse.FileType("r"), help='First Genbank file')
    parser.add_argument('gbk2', type=argparse.FileType("r"), help='Second Genbank file')

    args = parser.parse_args()

    cds1 = get_features_from_gbk(gbk_file=args.gbk1)
    cds2 = get_features_from_gbk(gbk_file=args.gbk2)
    print compare_feature_lists(list_a=cds1, list_b=cds2)
