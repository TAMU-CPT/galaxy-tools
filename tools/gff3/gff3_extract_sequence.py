#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_id(feature=None, parent_prefix=None):
    result = ""
    if parent_prefix is not None:
        result += parent_prefix + '|'
    if 'locus_tag' in feature.qualifiers:
        result += feature.qualifiers['locus_tag'][0]
    elif 'gene' in feature.qualifiers:
        result += feature.qualifiers['gene'][0]
    elif 'product' in feature.qualifiers:
        result += feature.qualifiers['product'][0]
    else:
        result += '%s_%s_%s' % (feature.location.start, feature.location.end,
                                feature.location.strand)
    return result


def feature_recurse(rec, feature, feature_filter):
    if type is None or feature.type == feature_filter:
        id = feature.id
        if len(id) == 0:
            id = get_id(feature)

        description = '[Location=%s]' % str(feature.location)
        yield [
            SeqRecord(
                feature.extract(rec).seq,
                id=id,
                description=description
            )
        ]
    if hasattr(feature, 'sub_features'):
        for subfeature in feature.sub_features:
            feature_recurse(rec, subfeature, feature_filter)


def main(fasta, gff3, feature_filter=None):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    for rec in GFF.parse(gff3, base_dict=seq_dict):
        for feat in rec.features:
            for hit in feature_recurse(rec, feat, feature_filter):
                yield hit

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export corresponding sequence in genome from GFF3', epilog="")
    parser.add_argument('fasta', type=file, help='Fasta Genome')
    parser.add_argument('gff3', help='GFF3 File')
    parser.add_argument('--feature_filter', default=None, help='Filter for specific feature types')
    args = parser.parse_args()

    for seq in main(**vars(args)):
        SeqIO.write(seq, sys.stdout, 'fasta')
