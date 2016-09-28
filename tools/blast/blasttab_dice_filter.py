#!/usr/bin/env python
import re
import sys
import copy
import argparse
from BCBio import GFF
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name='blasttab2gff3')

__doc__ = """
Blast TSV files, when transformed to GFF3, do not normally show gaps in the
blast hits. This tool aims to fill that "gap".
"""


def blasttsv2gff3(blasttsv, min_dice=50):
    columns = [
        'qseqid',      # 01 Query Seq-id (ID of your sequence)
        'sseqid',      # 02 Subject Seq-id (ID of the database hit)
        'pident',      # 03 Percentage of identical matches
        'length',      # 04 Alignment length
        'mismatch',    # 05 Number of mismatches
        'gapopen',     # 06 Number of gap openings
        'qstart',      # 07 Start of alignment in query
        'qend',        # 08 End of alignment in query
        'sstart',      # 09 Start of alignment in subject (database hit)
        'send',        # 10 End of alignment in subject (database hit)
        'evalue',      # 11 Expectation value (E-value)
        'bitscore',    # 12 Bit score
        'sallseqid',   # 13 All subject Seq-id(s), separated by a ';'
        'score',       # 14 Raw score
        'nident',      # 15 Number of identical matches
        'positive',    # 16 Number of positive-scoring matches
        'gaps',        # 17 Total number of gaps
        'ppos',        # 18 Percentage of positive-scoring matches
        'qframe',      # 19 Query frame
        'sframe',      # 20 Subject frame
        'qseq',        # 21 Aligned part of query sequence
        'sseq',        # 22 Aligned part of subject sequence
        'qlen',        # 23 Query sequence length
        'slen',        # 24 Subject sequence length
        'salltitles',  # 25 All subject title(s), separated by a '<>'
    ]

    print ' '.join(columns)
    for line in blasttsv:
        data = line.split('\t')
        dice = 2 * int(data[14]) / (float(data[22]) + float(data[23]))
        if dice >= min_dice:
            yield line


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Blast TSV to gapped GFF3')
    parser.add_argument('blasttsv', type=file, help='Blast TSV Output')
    parser.add_argument('--min_dice', type=int, help='Minimum dice score', default=0.5)
    args = parser.parse_args()

    for line in blasttsv2gff3(**vars(args)):
        print line
