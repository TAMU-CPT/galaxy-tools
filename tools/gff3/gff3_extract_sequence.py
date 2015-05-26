#!/usr/bin/env python
import argparse
from BCBio import GFF
from Bio import SeqIO

def feature_recurse(rec, feature, feature_filter):
    if type is None or feature.type == feature_filter:
        seq = str(feature.extract(rec).seq)
        print '>%s\n%s' % (feature.id, seq)
    if hasattr(feature, 'sub_features'):
        for subfeature in feature.sub_features:
            feature_recurse(rec, subfeature, feature_filter)


def main(fasta, gff3, feature_filter=None):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    codon_usage = {}

    for rec in GFF.parse(gff3, base_dict=seq_dict):
        for feat in rec.features:
            feature_recurse(rec, feat, feature_filter)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export corresponding sequence in genome from GFF3', epilog="")
    parser.add_argument('fasta', type=file, help='Fasta Genome')
    parser.add_argument('gff3', help='GFF3 File')
    parser.add_argument('--feature_filter', default=None, help='Filter for specific feature types')
    args = parser.parse_args()
    main(**vars(args))
