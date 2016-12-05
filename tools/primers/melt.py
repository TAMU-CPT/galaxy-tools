#!/usr/bin/env python
import argparse
from math import log10
import logging
from Bio import SeqIO
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name='melt')


def counts(sequence):
    return {
        'a': sequence.count('A'),
        't': sequence.count('T'),
        'c': sequence.count('C'),
        'g': sequence.count('G'),
    }


def nosalt(sequence, na=0.050):
    if len(sequence) < 14:
        def tm_func(a, c, t, g):
            return (a + t) * 2 + (g + c) * 4
    else:
        def tm_func(a, c, t, g):
            return 64.9 + 41 * (g + c - 16.4) / (a + c + t + g)

    return tm_func(**counts(sequence))


def salt(sequence, na=0.050):
    if len(sequence) < 14:
        def tm_func(a, c, t, g, na):
            return (a + t) * 2 + (g + c) * 4 - 16.6 * log10(0.050) + 16.6 * log10(na)
    else:
        def tm_func(a, c, t, g, na):
            return 100.5 + (41 * (g + c) / (a + t + g + c)) - (820 / (a + t + g + c)) + 16.6 * log10(na)

    kw = counts(sequence)
    kw['na'] = na

    return tm_func(**kw)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA melting point calculations')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='Fasta protein file')
    parser.add_argument('mode', choices=(
        'nosalt',
        'salt',
    ), default='nosalt')
    parser.add_argument('--na', type=float, help='Concentration of [Na+] in mM', default='50')
    args = parser.parse_args()

    tm_calculation = locals()[args.mode]

    print '#Sequence\tTm'
    for record in SeqIO.parse(args.fasta, "fasta"):
        print tm_calculation(str(record.seq).upper(), na=args.na / 1000)
