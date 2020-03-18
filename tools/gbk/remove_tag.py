#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
import sys
from Bio import SeqIO

import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def remove_qualifiers(genbank_files, feature_type=None, tag_type=None, tag_match=None):
    for genbank_file in genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            # Loop over features
            for feature in record.features:
                # If a feature_type is specified, use that.
                if feature_type is None or feature.type == feature_type:
                    # Loop over qualifiers (e.g. /locus_tag, /product)
                    good_qualifiers = {}
                    for key in feature.qualifiers.keys():
                        # If we haven't specified a tag_type (i.e. match any) or we
                        # have specified one
                        if tag_type is None or key == tag_type:
                            # Feature qualifiers can be multiply valued. They will
                            # always return a list of values. (e.g. /note="A", /note="B")
                            acceptable = []
                            for value in feature.qualifiers[key]:
                                # If tag_match is any, match all keys. If tag_match
                                # is in that qualifier, match that.
                                if tag_match is None or tag_match in value:
                                    pass
                                else:
                                    acceptable.append(value)
                            good_qualifiers[key] = acceptable
                        else:
                            good_qualifiers[key] = feature.qualifiers[key]

                    # Update the feature's qualifiers to "good" ones
                    feature.qualifiers = good_qualifiers

            # Print out the genbank file
            yield [record]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove specific features from a Genbank file"
    )
    parser.add_argument(
        "genbank_files", type=argparse.FileType("r"), nargs="+", help="Genbank files"
    )
    parser.add_argument("--feature_type", help="Feature type to remove")
    parser.add_argument("--tag_type", help="Tag type to remove")
    parser.add_argument("--tag_match", help="String in tag to match")

    args = parser.parse_args()
    for record in remove_qualifiers(**vars(args)):
        SeqIO.write(record, sys.stdout, "genbank")
