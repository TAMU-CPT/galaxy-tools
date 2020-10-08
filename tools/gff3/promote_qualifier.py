#!/usr/bin/env python
import argparse
import sys
import logging
from cpt_gffParser import gffParse, gffWrite
from gff3 import feature_lambda, feature_test_type

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def promote_qualifier(qualifier, parent, child, gff3):
    for record in gffParse(gff3):
        for parent_feature in feature_lambda(
            record.features, feature_test_type, {"type": parent}, subfeatures=True
        ):
            # for each feature of the parent type, get the first subfeature of the child type
            try:
                first_child = sorted(
                    list(
                        feature_lambda(
                            parent_feature.sub_features,
                            feature_test_type,
                            {"type": child},
                            subfeatures=False,
                        )
                    ),
                    key=lambda x: x.location.start
                    if parent_feature.strand > 0
                    else x.location.end,
                    reverse=False if parent_feature.strand > 0 else True,
                )[0]
            except IndexError:
                logging.warning("Child type %s not found under parent %s" % (child, parent_feature.qualifiers["ID"]))
                continue
            try:
                parent_feature.qualifiers[qualifier] = first_child.qualifiers[qualifier]
                logging.info(
                    "Promoted %s=%s in child %s to parent %s"
                    % (
                        qualifier,
                        first_child.qualifiers[qualifier],
                        first_child.qualifiers["ID"],
                        parent_feature.qualifiers["ID"],
                    )
                )
            except KeyError:
                logging.warning(
                    "Qualifier %s not found in child feature %s"
                    % (qualifier, first_child.qualifiers["ID"])
                )
        gffWrite([record], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Promote a child feature's qualifer to the parent feature's qualifier",
        epilog="",
    )
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 File")
    parser.add_argument(
        "parent",
        type=str,
        help="Feature type of the target parent feature (ex: gene, mrna, exon",
    )
    parser.add_argument(
        "child",
        type=str,
        help="Feature type of the target child feature (ex: mrna, exon, CDS",
    )
    parser.add_argument(
        "qualifier", help="Sepcific qualifier to promote (ex: Name, product, notes"
    )
    args = parser.parse_args()
    promote_qualifier(**vars(args))
