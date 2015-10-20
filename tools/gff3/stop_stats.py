#!/usr/bin/env python
import argparse
from BCBio import GFF
from Bio import SeqIO


def main(fasta, gff3):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    codon_usage = {}

    for rec in GFF.parse(gff3, base_dict=seq_dict):
        for feat in rec.features:
            seq = str(feat.extract(rec).seq)[-3:]
            try:
                codon_usage[seq] += 1
            except KeyError:
                codon_usage[seq] = 1

    names = {
        'TAG': 'Amber',
        'TAA': 'Ochre',
        'TGA': 'Opal',
    }

    # TODO: print all actg combinations? Or just ones that are there
    print '# Name\tCodon\tCount'
    for key in sorted(codon_usage):
        print '\t'.join((names.get(key.upper(), 'None'), key, str(codon_usage[key])))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarise stop codon usage', epilog="")
    parser.add_argument('fasta', type=file, help='Fasta Genome')
    parser.add_argument('gff3', help='GFF3 File')
    args = parser.parse_args()
    main(**vars(args))
