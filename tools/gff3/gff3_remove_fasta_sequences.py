#!/usr/bin/env python
"""
This program remove the fasta sequences from the end of a .gff3 file
"""

import argparse


def remove_fasta_seq(gff3, ogff3):
    i = 0
    stop = float("inf")
    # creates output file
    with open(ogff3, "w") as handle:
        # iterates ilne by line through input
        for line in gff3:
            i += 1
            # finds the '##FASTA' line, and line number becomes 'stop'
            if line.startswith("##FASTA"):
                stop = i
            # only writes lines to output if before stop
            if i < stop:
                handle.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Identify lysis gene candidates next to possible endolysin or holin genes",
        epilog="",
    )
    parser.add_argument("gff3", type=argparse.FileType("rw"), help="original gff3 file")
    parser.add_argument("--ogff3", type=str, default="output.gff3")
    args = parser.parse_args()

    remove_fasta_seq(args.gff3, args.ogff3)
