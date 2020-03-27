#!/usr/bin/env python
import sys
import logging
import argparse
from gff3 import feature_lambda, feature_test_type
from BCBio import GFF

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def safe_qualifiers(quals):
    unsafe_quals = ("ID", "Parent", "Name")
    new_quals = {}
    for (key, value) in quals.items():
        if key not in unsafe_quals:
            new_quals[key] = value

    return new_quals


def fix_ncbi(gff3):
    for rec in GFF.parse(gff3):
        for feature in feature_lambda(
            rec.features, feature_test_type, {"type": "gene"}, subfeatures=True
        ):
            CDSs = list(
                feature_lambda(
                    feature.sub_features,
                    feature_test_type,
                    {"type": "CDS"},
                    subfeatures=False,
                )
            )
            if len(CDSs) == 1:
                feature.qualifiers.update(safe_qualifiers(CDSs[0].qualifiers))

        GFF.write([rec], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="correct their terrible GFF data")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    fix_ncbi(**vars(args))
