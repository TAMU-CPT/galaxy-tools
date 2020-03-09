#!/usr/bin/env python
"""
This program remove the fasta sequences from the end of a .gff3 file
"""

import argparse


def remove_fasta_seq(gff3, ogff3):
    # iterates line by line through input
    for line in gff3:
        # writes lines to output until ##FASTA
        if line.startswith("##FASTA"):
            return
        else:
            ogff3.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Identify lysis gene candidates next to possible endolysin or holin genes",
        epilog="",
    )
    parser.add_argument("gff3", type=argparse.FileType("r"), help="original gff3 file")
    parser.add_argument("ogff3", type=argparse.FileType("w"), default="output.gff3")
    args = parser.parse_args()

    remove_fasta_seq(args.gff3, args.ogff3)
