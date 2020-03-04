#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import sys
import argparse
import copy
import logging
from BCBio import GFF
from Bio import SeqIO
from gff3 import feature_lambda, feature_test_true

logging.basicConfig(level=logging.INFO)


def extract_features(gff3_file):
    for rec in GFF.parse(gff3_file):
        for feat in feature_lambda(
            rec.features, feature_test_true, {}, subfeatures=False
        ):
            if feat.type == "remark":
                continue

            feat.qualifiers["color"] = ["255 0 0"]
            if feat.type == "Shine_Dalgarno_sequence":
                feat.type = "RBS"
                feat.qualifiers["color"] = ["180 0 0"]

            if feat.type not in ("CDS", "RBS", "gene", "terminator"):
                feat.type = "CDS"

            # Remove keys with '-'
            quals = copy.deepcopy(feat.qualifiers)
            for key in quals.keys():
                if "-" in key:
                    del quals[key]
            feat.qualifiers = quals

            yield feat


def merge_features(features=None, genbank_file=None):
    records = SeqIO.parse(genbank_file, "genbank")

    for record in records:
        for feature in features:
            record.features.append(feature)
        yield [record]


if __name__ == "__main__":
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description="Merge GFF3 data into a Genbank file")
    parser.add_argument(
        "genbank_file", type=argparse.FileType("r"), help="Genbank file"
    )
    parser.add_argument("gff3_file", type=argparse.FileType("r"), help="GFF3 Input")

    args = parser.parse_args()
    features = extract_features(args.gff3_file)
    for record in merge_features(features=features, genbank_file=args.genbank_file):
        SeqIO.write(record, sys.stdout, "genbank")
