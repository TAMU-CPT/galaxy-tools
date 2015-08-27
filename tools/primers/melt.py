#!/usr/bin/env python
import argparse
from Bio import SeqIO
import sys
from math import log10

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name='melt')

def counts(sequence):
    return {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'C': sequence.count('C'),
        'G': sequence.count('G'),
    }

def nosalt(sequence, na=0.050):
    if len(sequence) < 14:
        tm_func = "({A} + {T}) * 2 + ({G} + {C}) * 4"
    else:
        tm_func = "64.9 + 41 * ({G} + {C} - 16.4)/({A} + {T} + {C} + {G})"

    return eval(tm_func.format(**counts(sequence)))

def salt(sequence, na=0.050):
    if len(sequence) < 14:
        tm_func = '({A}+{T})*2 + ({G}+{C})*4 - 16.6*log10(0.050) + 16.6*log10({na})'
    else:
        tm_func = '100.5 + (41 * ({G}+{C})/({A}+{T}+{G}+{C})) - (820/({A}+{T}+{G}+{C})) + 16.6*log10({na})'

    return eval(tm_func.format(na=na, **counts(sequence)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA melting point calculations')
    parser.add_argument('fasta', type=file, help='Fasta protein file')
    parser.add_argument('mode', choices=(
        'nosalt',
        'salt',
    ), default='nosalt')
    parser.add_argument('--na', type=float, help='Concentration of [Na+] in mM', default='50')
    args = parser.parse_args()

    tm_calculation = locals()[args.mode]

    print '#Sequence\tTm'
    for record in SeqIO.parse(args.fasta, "fasta"):
        print tm_calculation(str(record.seq).upper(), na=args.na/1000)
