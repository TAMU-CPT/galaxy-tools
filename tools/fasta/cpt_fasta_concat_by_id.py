#!/usr/bin/env python
import argparse
from Bio import SeqIO
import logging
logging.basicConfig(level=logging.INFO)


def merge_fasta(left, right, rejects):
    matched = []
    left_side = list(SeqIO.parse(left, 'fasta'))
    for record in SeqIO.parse(right, 'fasta'):
        for left_record in left_side:
            if record.id == left_record.id:
                matched.append(left_record + record)
    import pprint; pprint.pprint(matched)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge two fasta files based on ID')
    parser.add_argument('left', type=file, help='Left hand side')
    parser.add_argument('right', help='Right hand side')
    parser.add_argument('--rejects', type=argparse.FileType('w'), help='Rejected matches', default='output.txt')

    args = parser.parse_args()
    merge_fasta(**vars(args))
