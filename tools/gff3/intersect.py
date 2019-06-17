#!/usr/bin/env python
import logging
import argparse
from intervaltree import IntervalTree, Interval
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def treeFeatures(features):
    for feat in features:
        # Interval(begin, end, data)
        yield Interval(int(feat.location.start), int(feat.location.end), feat.id)


def intersect(a, b):
    rec_a = list(GFF.parse(a))
    rec_b = list(GFF.parse(b))
    if len(rec_a) > 0 and len(rec_b) > 0:

        if len(rec_a) > 1 or len(rec_b) > 1:
            raise Exception("Cannot handle multiple GFF3 records in a file, yet")

        rec_a = rec_a[0]
        rec_b = rec_b[0]

        #builds interval tree from Interval objects of form (start, end, id) for each feature
        tree_a = IntervalTree(list(treeFeatures(rec_a.features)))
        tree_b = IntervalTree(list(treeFeatures(rec_b.features)))

        #Used to map ids back to features later
        rec_a_map = {f.id: f for f in rec_a.features}
        rec_b_map = {f.id: f for f in rec_b.features}

        rec_a_hits_in_b = []
        rec_b_hits_in_a = []

        for feature in rec_a.features:
            #Save each feature in rec_a that overlaps a feature in rec_b
            #hits = tree_b.find_range((int(feature.location.start), int(feature.location.end)))
            hits = tree_b[int(feature.location.start):int(feature.location.end)]
            #feature id is saved in interval result.data, use map to get full feature
            for hit in hits:
                rec_a_hits_in_b.append(rec_b_map[hit.data])

        for feature in rec_b.features:
            #Save each feature in rec_a that overlaps a feature in rec_b
            #hits = tree_a.find_range((int(feature.location.start), int(feature.location.end)))
            hits = tree_a[int(feature.location.start):int(feature.location.end)]
            #feature id is saved in interval result.data, use map to get full feature
            for hit in hits:
                rec_b_hits_in_a.append(rec_a_map[hit.data])

        #Remove duplicate features using sets
        rec_a.features = set(rec_a_hits_in_b)
        rec_b.features = set(rec_b_hits_in_a)

    else:
        #If one input is empty, output two empty result files
        rec_a = SeqRecord(Seq(''), "none")
        rec_b = SeqRecord(Seq(''), "none")
    return rec_a, rec_b


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('a', type=argparse.FileType("r"))
    parser.add_argument('b', type=argparse.FileType("r"))
    parser.add_argument('--oa', type=str, default='a_hits_in_b.gff')
    parser.add_argument('--ob', type=str, default='b_hits_in_a.gff')
    args = parser.parse_args()

    b, a = intersect(args.a, args.b)

    with open(args.oa, 'w') as handle:
        GFF.write([a], handle)

    with open(args.ob, 'w') as handle:
        GFF.write([b], handle)
