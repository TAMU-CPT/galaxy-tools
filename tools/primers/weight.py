#!/usr/bin/env python
import argparse
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name='melt')


def weight(seqobj):
    letters = list('ACTG')
    counts = {}
    mults = {
        'A': 313.21,
        'T': 304.2,
        'C': 289.18,
        'G': 329.21
    }

    seq = str(seqobj.seq).upper()
    for letter in letters:
        counts[letter] = seq.count(letter)

    mass = sum([counts[x] * mults[x] for x in letters]) - 61.96
    return [seqobj.id, mass]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA melting point calculations')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='Fasta protein file')
    args = parser.parse_args()

    print '#Sequence\tMass'
    for record in SeqIO.parse(args.fasta, "fasta"):
        print '\t'.join(map(str, weight(record)))
