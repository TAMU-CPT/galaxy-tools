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
    for idx, feature in enumerate(
        feature_lambda(
            rec.features, feature_test_type, {"type": "tRNA"}, subfeatures=True
        )
    ):
        fid = "tRNA-%03d" % (1 + idx)
        name = ["tRNA-" + feature.qualifiers["Codon"][0]]
        gene = SeqFeature(
            location=feature.location,
            type="gene",
            qualifiers={"ID": [fid + ".gene"], "source": ["aragorn"], "Name": name},
        )
        feature.qualifiers["Name"] = name
        # Below that we have an mRNA
        exon = SeqFeature(
            location=feature.location,
            type="exon",
            qualifiers={"source": ["aragorn"], "ID": ["%s.exon" % fid], "Name": name},
        )
        feature.qualifiers["ID"] = [fid]

        # gene -> trna -> exon
        feature.sub_features = [exon]
        gene.sub_features = [feature]
        yield gene


def gff_filter(gff3):
    for rec in GFF.parse(gff3):
        rec.features = sorted(list(fixed_feature(rec)), key=lambda x: x.location.start)
        rec.annotations = {}
        GFF.write([rec], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="add parent gene features to CDSs")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    gff_filter(**vars(args))
