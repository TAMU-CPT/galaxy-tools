#!/usr/bin/env python
import sys
import logging
import argparse
from Bio import SeqIO
from cpt_gffParser import gffParse, gffWrite
from gff3 import feature_lambda, feature_test_contains

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def gff_reopen(gff3, index=1, fasta=None, fasta_output=None):
    # Convert to zero-based
    index -= 1
    it = None
    if fasta:
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        it = gffParse(gff3, base_dict=seq_dict)
    else:
        it = gffParse(gff3)

    for rec in it:
        # Reopen
        if (
            len(
                list(
                    feature_lambda(
                        rec.features,
                        feature_test_contains,
                        {"index": index},
                        subfeatures=False,
                    )
                )
            )
            > 0
        ):
            log.warn(
                "WARNING: Index chosen is in the middle of a feature. This feature will disappear from the output"
            )
        # TODO: This call removes metadata!
        rec = rec[index:] + rec[0:index]
        rec.features = sorted(rec.features, key=lambda x: x.location.start)

        if fasta:
            if len(rec.seq) == rec.seq.count("?"):
                log.error(
                    "ERROR: You have provided a fasta file but the sequence ID in the fasta file DID NOT MATCH THE GFF. THIS IS BAD."
                )

        yield rec


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reopen a set of GFF3 annotations")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    parser.add_argument(
        "--fasta", type=argparse.FileType("r"), help="Optional fasta file"
    )
    parser.add_argument(
        "--fasta_output",
        type=argparse.FileType("w"),
        help="Optional fasta file output",
        default="reopened.fasta",
    )
    parser.add_argument("index", type=int, help="Index to reopen genome at")
    args = parser.parse_args()

    for rec in gff_reopen(**vars(args)):
        gffWrite([rec], sys.stdout)
        if args.fasta:
            SeqIO.write([rec], args.fasta_output, "fasta")
