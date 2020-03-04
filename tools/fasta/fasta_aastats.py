#!/usr/bin/env python
import logging
import argparse
import string
from Bio import SeqIO

logging.basicConfig(level=logging.INFO)


def aa_stats(fasta_file, **kwargs):
    records = list(SeqIO.parse(fasta_file, "fasta"))

    tn_table_keys = list(string.ascii_uppercase) + ["*"]

    header = ["#ID", "Length"] + tn_table_keys
    yield header

    for record in records:
        aa_counts = {}

        for nt in list(str(record.seq)):
            try:
                aa_counts[nt] += 1
            except:
                aa_counts[nt] = 1

        row = [record.id, len(record.seq)]
        numbers = []
        for nt in tn_table_keys:
            if nt in aa_counts:
                numbers.append(aa_counts[nt])
            else:
                numbers.append(0)

        numbers = [float(x) / sum(numbers) for x in numbers]
        yield row + numbers


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate AA frequencies in sequences"
    )
    parser.add_argument("fasta_file", type=argparse.FileType("r"), help="Fasta file")
    parser.add_argument("--version", action="version", version="0.1")
    args = parser.parse_args()

    for row in aa_stats(**vars(args)):
        print("\t".join(map(str, row)))
