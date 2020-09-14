#!/usr/bin/env python
import sys
import logging
import argparse
from cpt_gffParser import gffParse, gffWrite
from gff3 import feature_lambda

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def test_true(feature, **kwargs):
    return True


def gff_filter(gff3):
    for rec in gffParse(gff3):
        for feature in feature_lambda(rec.features, test_true, {}, subfeatures=True):
            if feature.type == "exon" and len(feature) < 20:
                feature.type = "Shine_Dalgarno_sequence"

        rec.annotations = {}
        gffWrite([rec], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="original annotation set")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    gff_filter(**vars(args))
