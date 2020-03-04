#!/usr/bin/env python
import sys
import logging
import argparse
from Bio import SeqIO
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def find_introns(gff3, fasta):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    for rec in GFF.parse(gff3, base_dict=seq_dict):
        genes = list(
            feature_lambda(
                rec.features, feature_test_type, {"type": "gene"}, subfeatures=True
            )
        )
        for gene in genes:
            cdss = sorted(
                list(
                    feature_lambda(
                        gene.sub_features,
                        feature_test_type,
                        {"type": "CDS"},
                        subfeatures=False,
                    )
                ),
                key=lambda x: x.location.start,
            )
            if len(cdss) > 1:
                intron = ""
                for i in range(
                    len(cdss) - 1
                ):  # find pairs of cdss with introns in between
                    intron_start = cdss[i].location.end
                    intron_end = cdss[i + 1].location.start
                    intron += rec[intron_start:intron_end].seq
                sys.stdout.write(">" + rec.id + "\n")
                sys.stdout.write(intron + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="extract introns from gene features that have more than one CDS"
    )
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    parser.add_argument("fasta", type=argparse.FileType("r"), help="fasta file")
    args = parser.parse_args()
    find_introns(**vars(args))
