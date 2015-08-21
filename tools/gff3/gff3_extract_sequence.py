#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gff3 import feature_lambda, feature_test_type, get_id


def main(fasta, gff3, feature_filter=None):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    for rec in GFF.parse(gff3, base_dict=seq_dict):
        for feat in feature_lambda(
            rec.features,
            feature_test_type,
            {'type': feature_filter},
            subfeatures = False
        ):
            id = feat.id
            if len(id) == 0:
                id = get_id(feat)

            description = '[Location=%s]' % str(feat.location)
            yield [
                SeqRecord(
                    feat.extract(rec).seq,
                    id=id,
                    description=description
                )
            ]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export corresponding sequence in genome from GFF3', epilog="")
    parser.add_argument('fasta', type=file, help='Fasta Genome')
    parser.add_argument('gff3', help='GFF3 File')
    parser.add_argument('--feature_filter', default=None, help='Filter for specific feature types')
    args = parser.parse_args()

    for seq in main(**vars(args)):
        SeqIO.write(seq, sys.stdout, 'fasta')
