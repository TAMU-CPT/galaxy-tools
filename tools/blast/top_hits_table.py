#!/usr/bin/env python
import sys
import os
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type, fsort, get_id
from relatedness import parse_blast, with_dice


def important_only(blast):
    for data in blast:
        yield [
            data[0],  # 01 Query Seq-id (ID of your sequence)
            # 02 Subject Seq-id (ID of the database hit)
            # 03 Percentage of identical matches
            # 04 Alignment length
            # 05 Number of mismatches
            # 06 Number of gap openings
            # 07 Start of alignment in query
            # 08 End of alignment in query
            # 09 Start of alignment in subject (database hit)
            # 10 End of alignment in subject (database hit)
            float(data[10]),  # 11 Expectation value (E-value)
            # 12 Bit score
            # 13 All subject Seq-id(s), separated by a ';'
            zip(
                data[12].split(';'),
                data[24].split('<>')
            ),
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
            data[25]  # 26 dice
        ]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gff3', type=argparse.FileType("r"), help='Query Genome Features')
    parser.add_argument('blasttsv', type=argparse.FileType("r"), nargs='+', help='Blast TSV Output')
    args = parser.parse_args()

    top_hits = {}
    blast_names = []
    for fh in args.blasttsv:
        fn = os.path.basename(fh.name)
        blast_names.append(fn)
        data = important_only(with_dice(parse_blast(fh)))
        for (qseq, evalue, sseq, dice) in data:
            if qseq not in top_hits:
                top_hits[qseq] = {}

            if fn not in top_hits[qseq]:
                top_hits[qseq][fn] = (evalue, sseq, dice)

            if evalue <= top_hits[qseq][fn][0] and dice >= top_hits[qseq][fn][2]:
                top_hits[qseq][fn] = (evalue, sseq, dice)


    sys.stdout.write('# Query Feature\tLocation\t')
    sys.stdout.write('\t'.join(['%s\tevalue\tdice' % x for x in blast_names]))
    sys.stdout.write('\n')
    for rec in GFF.parse(args.gff3):
        for feat in fsort(feature_lambda(
            rec.features,
            feature_test_type,
            {'types': 'CDS'},
            subfeatures=False,
        )):
            sys.stdout.write(feat._parent._parent.qualifiers['Name'][0])
            sys.stdout.write('\t')
            sys.stdout.write(str(feat.location))

            for db in blast_names:
                fid = get_id(feat)
                if fid in top_hits:
                    if fn in top_hits[fid]:
                        sys.stdout.write('\t')
                        sys.stdout.write(';'.join(['%s %s' % (x, y) for (x, y) in top_hits[fid][fn][1]]))
                        sys.stdout.write('\t')
                        sys.stdout.write(str(top_hits[fid][fn][0]))
                        sys.stdout.write('\t')
                        sys.stdout.write(str(top_hits[fid][fn][2]))
                    else:
                        sys.stdout.write('\tNone')
                        sys.stdout.write('\tNone')
                        sys.stdout.write('\tNone')
                else:
                    sys.stdout.write('\tNone')
                    sys.stdout.write('\tNone')
                    sys.stdout.write('\tNone')

            sys.stdout.write('\n')
