#!/usr/bin/env python
# import sys
import argparse
from BCBio import GFF
from Bio import SeqIO

def taper_list(gff3, fasta):
    """ delete records that already have identified transmembrane domains """
    records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    for rec in GFF.parse(gff3):
        del records[rec.id]

    find_tmembrane(records)


def find_tmembrane(records):
    """ identify transmembrane domains based on the following rules:
            (1) domain must either have 16 consecutive neutrally-charged residues or
            (2) if (1) is not met, allow for lysine (K) residues at locations n to n+2 and/or n+13 to n+15
    """

    return



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='find phage transmembrane domains')
    parser.add_argument('gff3', type=file, help='GFF3 output of TMHMM')
    parser.add_argument('fasta', type=file, help='fasta file of protein(s)')
    args = parser.parse_args()
    taper_list(**vars(args))
