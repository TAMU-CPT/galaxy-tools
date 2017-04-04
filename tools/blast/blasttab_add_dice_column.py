#!/usr/bin/env python
import sys
import argparse


def blasttsv2gff3(blasttsv):
    # 01 Query Seq-id (ID of your sequence)
    # 02 Subject Seq-id (ID of the database hit)
    # 03 Percentage of identical matches
    # 04 Alignment length
    # 05 Number of mismatches
    # 06 Number of gap openings
    # 07 Start of alignment in query
    # 08 End of alignment in query
    # 09 Start of alignment in subject (database hit)
    # 10 End of alignment in subject (database hit)
    # 11 Expectation value (E-value)
    # 12 Bit score
    # 13 All subject Seq-id(s), separated by a ';'
    # 14 Raw score
    # 15 Number of identical matches
    # 16 Number of positive-scoring matches
    # 17 Total number of gaps
    # 18 Percentage of positive-scoring matches
    # 19 Query frame
    # 20 Subject frame
    # 21 Aligned part of query sequence
    # 22 Aligned part of subject sequence
    # 23 Query sequence length
    # 24 Subject sequence length
    # 25 All subject title(s), separated by a '<>'

    for line in blasttsv:
        data = line.split('\t')
        dice = 2 * int(data[14]) / (float(data[22]) + float(data[23]))
        yield line.strip('\n') + '\t' + str(dice)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Blast TSV to gapped GFF3')
    parser.add_argument('blasttsv', type=argparse.FileType("r"), help='Blast TSV Output')
    args = parser.parse_args()

    for line in blasttsv2gff3(**vars(args)):
        sys.stdout.write(line)
        sys.stdout.write('\n')
