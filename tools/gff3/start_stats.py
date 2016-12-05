#!/usr/bin/env python
import argparse
from BCBio import GFF
from Bio import SeqIO
from gff3 import feature_lambda, feature_test_type


def main(fasta, gff3):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    codon_usage = {}

    for rec in GFF.parse(gff3, base_dict=seq_dict):
        for feat in feature_lambda(rec.features, feature_test_type, {'type': 'CDS'}, subfeatures=True):
            seq = str(feat.extract(rec).seq)[0:3]
            try:
                codon_usage[seq] += 1
            except KeyError:
                codon_usage[seq] = 1

    # TODO: print all actg combinations? Or just ones that are there
    print '# Codon\tCount'
    for key in sorted(codon_usage):
        print '\t'.join((key, str(codon_usage[key])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarise start codon usage', epilog="")
    parser.add_argument('fasta', type=argparse.FileType("r"), help='Fasta Genome')
    parser.add_argument('gff3', help='GFF3 File')
    args = parser.parse_args()
    main(**vars(args))
