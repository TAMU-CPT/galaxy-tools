#!/usr/bin/env python
import sys
import copy
import logging
from interval_tree import IntervalTree
logging.basicConfig(level=logging.INFO)
import argparse
from gff3 import feature_lambda, feature_test_type
from BCBio import GFF
from Bio.SeqFeature import FeatureLocation
log = logging.getLogger(__name__)


def __get_features(child, interpro=False):
    child_features = {}
    for rec in GFF.parse(child):
        log.info("Parsing %s", rec.id)
        for feature in rec.features:
            parent_feature_id = rec.id
            if interpro:
                if feature.type == 'polypeptide':
                    continue
                if '_' in parent_feature_id:
                    parent_feature_id = parent_feature_id[parent_feature_id.index('_') + 1:]

            try:
                child_features[parent_feature_id].append(feature)
            except KeyError:
                child_features[parent_feature_id] = [feature]
    return child_features


def __update_feature_location(feature, parent, protein2dna):
    start = feature.location.start
    end = feature.location.end
    if protein2dna:
        start *= 3
        end *= 3

    if parent.location.strand >= 0:
        ns = parent.location.start + start
        ne = parent.location.start + end
        st = +1
    else:
        ns = parent.location.end - end
        ne = parent.location.end - start
        st = -1

    # Don't let start/stops be less than zero. It's technically valid for them
    # to be (at least in the model I'm working with) but it causes numerous
    # issues.
    #
    # Instead, we'll replace with %3 to try and keep it in the same reading
    # frame that it should be in.
    if ns < 0:
        ns %= 3
    if ne < 0:
        ne %= 3

    feature.location = FeatureLocation(ns, ne, strand=st)

    if hasattr(feature, 'sub_features'):
        for subfeature in feature.sub_features:
            __update_feature_location(subfeature, parent, protein2dna)

def treeFeatures(features, strand=0):
    for feat in features:
        if feat.location.strand == strand:
            yield (int(feat.location.start), int(feat.location.end), feat.id)


def neighbours(a, b, within=1000, mode='unordered'):
    rec_a = list(GFF.parse(a))
    rec_b = list(GFF.parse(b))
    if len(rec_a) > 1 or len(rec_b) > 1:
        raise Exception("Cannot handle multiple GFF3 records in a file, yet")

    rec_a = rec_a[0]
    rec_b = rec_b[0]

    tree_f = IntervalTree(list(treeFeatures(rec_a.features, strand=1)), 1, len(rec_a))
    tree_r = IntervalTree(list(treeFeatures(rec_a.features, strand=-1)), 1, len(rec_a))

    for feature in rec_b.features:
        start = feature.location.start
        end = feature.location.end
        if feature.location.strand > 0:
            start -= within
            if mode != 'ordered':
                end += within

            hits = tree_f.find_range((start, end))
            if len(hits) == 0: continue

            print start, end, feature.location.strand, hits, feature.id
        else:
            end += within
            if mode != 'ordered':
                start -= within

            hits = tree_r.find_range((start, end))
            if len(hits) == 0: continue

            print start, end, feature.location.strand, hits, feature.id


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('a', type=file)
    parser.add_argument('b', type=file)
    parser.add_argument('--within', type=int, default=1000)
    parser.add_argument('--mode', type=str, choices=('ordered', 'unordered'), default='unordered')
    args = parser.parse_args()

    neighbours(**vars(args))
