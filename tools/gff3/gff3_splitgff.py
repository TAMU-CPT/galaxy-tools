#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from cpt_gffParser import gffParse, gffWrite

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sample script to add an attribute to a feature via web services"
    )
    parser.add_argument("data", type=argparse.FileType("r"), help="GFF3 File")
    parser.add_argument(
        "--gff",
        type=argparse.FileType("w"),
        help="Output Annotations",
        default="data.gff3",
    )
    parser.add_argument(
        "--fasta",
        type=argparse.FileType("w"),
        help="Output Sequence",
        default="data.fa",
    )
    args = parser.parse_args()

    for record in gffParse(args.data):
        gffWrite([record], args.gff)
        record.description = ""
        if not isinstance(record.seq, str):
          record.seq = str(record.seq)
        SeqIO.write([record], args.fasta, "fasta")
        sys.exit()
