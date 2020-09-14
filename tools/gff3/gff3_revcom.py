#!/usr/bin/env python
import sys
import logging
import argparse
from Bio import SeqIO
from cpt_gffParser import gffParse, gffWrite

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def gff_filter(gff3, fasta=None, fasta_output=None):
    if fasta:
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        it = gffParse(gff3, base_dict=seq_dict)
    else:
        it = gffParse(gff3)

    for rec in it:
        rec.annotations = {}
        rec = rec.reverse_complement(
            id=True,
            name=True,
            description=True,
            features=True,
            annotations=True,
            letter_annotations=True,
            dbxrefs=True,
        )
        yield rec


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="reverse and complement as set of GFF3 annotations"
    )
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
    args = parser.parse_args()

    for rec in gff_filter(**vars(args)):
        gffWrite([rec], sys.stdout)
        if args.fasta:
            SeqIO.write([rec], args.fasta_output, "fasta")
