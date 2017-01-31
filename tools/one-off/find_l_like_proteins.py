#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO


def find_hydrophobic_seq(seq):
    hydrophobic_residues = 'A G V L I P F H S T Y C M W Q N'.split()  # allowed hydrophobic amino acids

    count = 0  # count for number of hydrophobic residues in a row
    start = 0
    end = 0
    domain = ''
    for num, s in enumerate(seq):
        if s in hydrophobic_residues:
            if not count:
                start = num
            count += 1
        else:
            print '***'
            print count
            end = num
            domain = seq[start:end]
            print {'start': start, 'end': end, 'domain': domain}
            print '***'
            count = 0
    sys.exit()


def find_l_like_proteins(fasta):
    """ Returns proteins that have:
            >= 10 hydrophobic aa's in a row (AGVLIPFHSTYCMWQN) TODO: QN pref at the ends
            'LS' right after that domain
            minimum net charge of +2 for N terminus
    """
    records = list(SeqIO.parse(fasta, "fasta"))

    for record in records:
        find_hydrophobic_seq(record.seq)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Finds L-like proteins')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='Fasta file')

    args = parser.parse_args()
    find_l_like_proteins(args.fasta)
