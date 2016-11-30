#!/usr/bin/env python
import logging
import random
import argparse
import StringIO
from Bio import SeqIO, Seq
logging.basicConfig(level=logging.INFO)


def disambiguate(fasta_file, seed=42, tbl_out=None):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    output = StringIO.StringIO()

    if seed != 0:
        random.seed(seed)

    replace = {
        'U': 'T',
        'R': 'AG',
        'Y': 'CT',
        'M': 'CA',
        'K': 'TG',
        'W': 'TA',
        'S': 'CG',
        'B': 'CTG',
        'D': 'ATG',
        'H': 'ATC',
        'V': 'ACG',
        'N': 'ACTG',
        '*': 'ACTG',
    }
    delta_tbl = [('# Pos', 'Orig', 'New')]

    for record in records:
        replacement_seq = ""
        for i, char in enumerate(str(record.seq).upper()):
            if char in replace:
                newchar = random.choice(replace[char])
                delta_tbl.append((i, char, newchar))
                replacement_seq += newchar
            elif char not in replace and char not in 'ACTG':
                delta_tbl.append((i, char, 'Unknown'))
                replacement_seq += char
            else:
                replacement_seq += char

        record.seq = Seq.Seq(replacement_seq)
        SeqIO.write(record, output, "fasta")

    with open(tbl_out, 'w') as handle:
        for row in delta_tbl:
            handle.write('\t'.join(map(str, row)) + "\n")

    print output.getvalue()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Replace ambiguous bases')
    parser.add_argument('fasta_file', type=file, help='Fasta file')
    parser.add_argument('--tbl_out', type=str, help='Table output', default='ambiguities.tsv')
    parser.add_argument('--seed', type=int, help='Seed for reproducible analysis', default=42)

    args = parser.parse_args()
    disambiguate(**vars(args))
