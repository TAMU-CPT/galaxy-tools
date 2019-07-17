#!/usr/bin/env python
import sys
import logging
import argparse
from gff3 import feature_lambda, feature_test_qual_value
from BCBio import GFF
from Bio.SeqFeature import FeatureLocation
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def __get_features(child, interpro=False):
    child_features = {}
    for rec in GFF.parse(child):
        log.info("Parsing %s", rec.id)
        # Only top level
        for feature in rec.features:
            # Get the record id as parent_feature_id (since this is how it will be during remapping)
            parent_feature_id = rec.id
            # If it's an interpro specific gff3 file
            if interpro:
                # Then we ignore polypeptide features as they're useless
                if feature.type == 'polypeptide':
                    continue
               
            try:
                child_features[parent_feature_id].append(feature)
            except KeyError:
                child_features[parent_feature_id] = [feature]
            # Keep a list of feature objects keyed by parent record id
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
   
    # Don't let start/stops be less than zero.
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


def rebase(parent, child, interpro=False, protein2dna=False, map_by='ID'):
    # get all of the features we will be re-mapping in a dictionary, keyed by parent feature ID
    child_features = __get_features(child, interpro=interpro)

    for rec in GFF.parse(parent):
        replacement_features = []
        # Horrifically slow I believe
        for feature in feature_lambda(
                rec.features,
                # Filter features in the parent genome by those that are
                # "interesting", i.e. have results in child_features array.
                # Probably an unnecessary optimisation.
                feature_test_qual_value,
                {
                    'qualifier': map_by,
                    'attribute_list': child_features.keys(),
                },
                subfeatures=False):

            # Features which will be re-mapped
            to_remap = child_features[feature.id]
            
            fixed_features = []
            for x in to_remap:
                # Then update the location of the actual feature
                __update_feature_location(x, feature, protein2dna)

                if interpro:
                    for y in ('status', 'Target'):
                        try:
                            del x.qualifiers[y]
                        except:
                            pass

                fixed_features.append(x)
            replacement_features.extend(fixed_features)
        # We do this so we don't include the original set of features that we
        # were rebasing against in our result.
        rec.features = replacement_features
        rec.annotations = {}
        GFF.write([rec], sys.stdout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('parent', type=argparse.FileType("r"), help='Parent GFF3 annotations')
    parser.add_argument('child', help='Child GFF3 annotations to rebase against parent')
    parser.add_argument('--interpro', action='store_true',
                        help='Interpro specific modifications')
    parser.add_argument('--protein2dna', action='store_true',
                        help='Map protein translated results to original DNA data')
    parser.add_argument('--map_by', help='Map by key', default='ID')
    args = parser.parse_args()
    rebase(**vars(args))
