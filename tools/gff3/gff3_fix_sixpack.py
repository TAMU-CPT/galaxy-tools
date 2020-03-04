#!/usr/bin/env python
import sys
import logging
import argparse
from BCBio import GFF
from Bio.SeqFeature import SeqFeature
from gff3 import feature_lambda, feature_test_type

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def fixed_feature(rec):
    # Get all gene features to remove the mRNAs from
    for feature in feature_lambda(
        rec.features, feature_test_type, {"type": "gene"}, subfeatures=True
    ):
        gene = feature
        sub_features = []
        # Filter out mRNA subfeatures, save other ones to new gene object.
        for sf in feature_lambda(
            feature.sub_features,
            feature_test_type,
            {"type": "mRNA"},
            subfeatures=True,
            invert=True,
        ):
            sf.qualifiers["Parent"] = gene.qualifiers["ID"]
            sub_features.append(sf)
        # override original subfeatures with our filtered list
        gene.sub_features = sub_features
        yield gene


def gff_filter(gff3):
    for rec in GFF.parse(gff3):
        rec.features = sorted(list(fixed_feature(rec)), key=lambda x: x.location.start)
        rec.annotations = {}
        GFF.write([rec], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix gene model from naive ORF caller")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    gff_filter(**vars(args))
