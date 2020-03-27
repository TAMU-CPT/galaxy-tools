#!/usr/bin/env python
import argparse
from Bio import SeqIO

import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def fasta_join(left, right):
    right_record_list = list(SeqIO.parse(right, "fasta"))
    right_records = {record.id: record for record in right_record_list}

    # has left right missing
    right_missing = open("right_missing.fa", "w")
    both_left_seq = open("both_left.fa", "w")
    both_right_seq = open("both_right.fa", "w")
    # has right missing left
    left_missing = open("left_missing.fa", "w")

    for record in SeqIO.parse(left, "fasta"):
        if record.id in right_records:
            SeqIO.write(record, both_left_seq, "fasta")
            SeqIO.write(right_records[record.id], both_right_seq, "fasta")
            del right_records[record.id]
        else:
            SeqIO.write(record, right_missing, "fasta")

    for record in right_record_list:
        if record.id in right_records:
            SeqIO.write(record, left_missing, "fasta")

    right_missing.close()
    both_left_seq.close()
    both_right_seq.close()
    left_missing.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "left", type=argparse.FileType("r"), help="For all sequences in LEFT"
    )
    parser.add_argument(
        "right", type=argparse.FileType("r"), help="Find matches in RIGHT"
    )

    args = parser.parse_args()
    fasta_join(**vars(args))
