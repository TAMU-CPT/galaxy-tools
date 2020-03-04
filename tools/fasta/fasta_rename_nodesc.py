#!/usr/bin/env python
import sys
import logging
from Bio import SeqIO
import argparse

logging.basicConfig(level=logging.INFO)


def rename_fasta_sequences(fasta_file, new_name):
    for idx, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        if idx > 1:
            raise Exception("Too many sequences")

        record.id = new_name.strip()
        record.description = ""

        SeqIO.write([record], sys.stdout, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rename fasta sequences and fully removes previous name"
    )
    parser.add_argument(
        "fasta_file",
        metavar="N",
        type=argparse.FileType("r"),
        nargs="?",
        help="fasta file",
    )
    parser.add_argument("new_name", nargs="?", help="New name for the fasta sequence")
    args = parser.parse_args()

    rename_fasta_sequences(**vars(args))
