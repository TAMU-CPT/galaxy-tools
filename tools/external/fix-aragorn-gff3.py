#!/usr/bin/env python
import sys
import logging
import argparse
from cpt_gffParser import gffParse, gffWrite, gffSeqFeature
from Bio.SeqFeature import SeqFeature
from gff3 import feature_lambda, feature_test_type

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def fixed_feature(rec):
    for idx, feature in enumerate(
        feature_lambda(
            rec.features, feature_test_type, {"types": ["tRNA", "tmRNA"]}, subfeatures=True
        )
    ):
        
        fid = "%s-%03d" % (feature.type, 1 + idx)
        try:
            name = [feature.type + "-" + feature.qualifiers["Codon"][0]]
        except KeyError:
            name = [feature.qualifiers['product'][0]]
        try:
          origSource = feature.qualifiers["source"][0]
        except:
          origSource = "."
        gene = gffSeqFeature(
            location=feature.location,
            type="gene",
            qualifiers={"ID": [fid + ".gene"], "source": [origSource], "Name": name},
        )
        feature.qualifiers["Name"] = name
        # Below that we have an mRNA
        exon = gffSeqFeature(
            location=feature.location,
            type="exon",
            qualifiers={"source": [origSource], "ID": ["%s.exon" % fid], "Name": name},
        )
        feature.qualifiers["ID"] = [fid]
        exon.qualifiers["Parent"] = [fid]
        feature.qualifiers["Parent"] = [fid + ".gene"]
        # gene -> trna -> exon
        feature.sub_features = [exon]
        gene.sub_features = [feature]
        yield gene


def gff_filter(gff3):
    found_gff = False
    for rec in gffParse(gff3):
        found_gff = True
        rec.features = sorted(list(fixed_feature(rec)), key=lambda x: x.location.start)
        rec.annotations = {}
        gffWrite([rec], sys.stdout)
    if not found_gff:
        print("##gff-version 3")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="add parent gene features to CDSs")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    gff_filter(**vars(args))
