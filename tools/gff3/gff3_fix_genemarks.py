#!/usr/bin/env python
import sys
import copy
import logging
import argparse
from cpt_gffParser import gffParse, gffWrite

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def fix_genemarks(gff3):
    for rec in gffParse(gff3):
        new_features = []
        for feature in rec.features:
            gene_parent = copy.deepcopy(feature)
            gene_parent.type = "gene"
            feature.qualifiers["ID"][0] = feature.qualifiers["ID"][0].replace(
                "gene", "cds"
            )
            gene_parent.sub_features = [feature]
            del gene_parent.qualifiers["score"]
            del gene_parent.qualifiers["phase"]
            new_features.append(gene_parent)
        rec.features = new_features
        gffWrite([rec], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="correct their terrible GFF data")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    fix_genemarks(**vars(args))
