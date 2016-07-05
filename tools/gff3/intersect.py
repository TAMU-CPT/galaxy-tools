#!/usr/bin/env python
import logging
import argparse
from interval_tree import IntervalTree
from BCBio import GFF
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def treeFeatures(features):
    for feat in features:
        yield (int(feat.location.start), int(feat.location.end), feat.id)


def intersect(a, b):
    rec_a = list(GFF.parse(a))
    rec_b = list(GFF.parse(b))
    if len(rec_a) > 1 or len(rec_b) > 1:
        raise Exception("Cannot handle multiple GFF3 records in a file, yet")

    rec_a = rec_a[0]
    rec_b = rec_b[0]

    tree_a = IntervalTree(list(treeFeatures(rec_a.features)), 1, len(rec_a))
    tree_b = IntervalTree(list(treeFeatures(rec_b.features)), 1, len(rec_b))

    rec_a_map = {f.id: f for f in rec_a.features}
    rec_b_map = {f.id: f for f in rec_b.features}

    rec_a_hits_in_b = []
    rec_b_hits_in_a = []

    for feature in rec_a.features:
        hits = tree_b.find_range((int(feature.location.start), int(feature.location.end)))
        for hit in hits:
            rec_a_hits_in_b.append(rec_b_map[hit])

    for feature in rec_b.features:
        hits = tree_a.find_range((int(feature.location.start), int(feature.location.end)))
        for hit in hits:
            rec_b_hits_in_a.append(rec_a_map[hit])

    rec_a.features = rec_a_hits_in_b
    rec_b.features = rec_b_hits_in_a
    return rec_a, rec_b


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('a', type=file)
    parser.add_argument('b', type=file)
    parser.add_argument('--oa', type=str, default='a_hits_in_b.gff')
    parser.add_argument('--ob', type=str, default='b_hits_in_a.gff')
    args = parser.parse_args()

    b, a = intersect(args.a, args.b)

    with open(args.oa, 'w') as handle:
        GFF.write([a], handle)

    with open(args.ob, 'w') as handle:
        GFF.write([b], handle)
