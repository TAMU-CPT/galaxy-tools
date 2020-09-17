#!/usr/bin/env python
import sys
import logging
import argparse
from cpt_gffParser import gffParse, gffWrite, gffSeqFeature
#from Bio.SeqFeature import SeqFeature
from gff3 import feature_lambda, feature_test_type

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def fixed_feature(rec):
    for feature in feature_lambda(
        rec.features, feature_test_type, {"type": "CDS"}, subfeatures=True
    ):
        import random

        fid = feature.qualifiers["ID"][0] + "_" + str(random.random())
        gene = gffSeqFeature(
            location=feature.location,
            type="gene",
            qualifiers={"ID": [fid], "source": ["cpt.fixModel"]},
        )
        # Below that we have an mRNA
        mRNA = gffSeqFeature(
            location=feature.location,
            type="mRNA",
            qualifiers={"source": ["cpt.fixModel"], "ID": ["%s.mRNA" % fid]},
        )
        feature.qualifiers["ID"] = [fid + ".CDS"]

        mRNA.sub_features = [feature]
        gene.sub_features = [mRNA]
        yield gene


def gff_filter(gff3):
    for rec in gffParse(gff3):
        rec.features = sorted(list(fixed_feature(rec)), key=lambda x: x.location.start)
        rec.annotations = {}
        gffWrite([rec], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="add parent gene features to CDSs")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    gff_filter(**vars(args))
