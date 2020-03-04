#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
import sys
import logging
from Bio import SeqIO

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def findById(genbank_files, feature_type=None, tag_type=None, tag_match=None):
    for genbank_file in genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            # Loop over features
            good_features = []
            for feature in record.features:
                # If a feature_type is specified, use that.
                if feature_type is not None and feature.type != feature_type:
                    continue

                for key in feature.qualifiers.keys():
                    # If we haven't specified a tag_type (i.e. match any) or we
                    # have specified one
                    if tag_type is not None and tag_type != key:
                        continue

                    for value in feature.qualifiers.get(key, []):
                        # If tag_match is any, match all keys. If tag_match
                        # is in that qualifier, match that.
                        if tag_match is not None and tag_match not in value:
                            continue

                        good_features.append(feature)

            record.features = good_features
            # Print out the genbank file
            yield [record]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove specific features from a Genbank file"
    )
    parser.add_argument(
        "genbank_files", type=argparse.FileType("r"), nargs="+", help="Genbank files"
    )
    parser.add_argument("--feature_type", help="Feature type")
    parser.add_argument("--tag_type", help="Tag type")
    parser.add_argument("--tag_match", help="String in tag")

    args = parser.parse_args()
    for record in findById(**vars(args)):
        SeqIO.write(record, sys.stdout, "genbank")
