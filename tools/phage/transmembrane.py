#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
from Bio import SeqIO

num_caught = 0

def taper_list(gff3, fasta):
    """ delete records that already have identified transmembrane domains """
    records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    for rec in GFF.parse(gff3):
        del records[rec.id]

    find_tmembrane(records)

def hydrophobicity(beg, mid, end):
    """ analyzes hydrophobicity of a sequence """
    # neutral amino acids
    aas = ['F', 'I', 'W', 'L', 'V', 'M', 'Y', 'C', 'A', 'T', 'G', 'S']

    for i in mid:
        if i not in aas:
            return False
        else:
            for j in beg+end:
                if j not in aas and j != 'K':
                    return False
    return True

def print_seq(plot, seq):
    print seq
    annotate = ""
    for i in range(len(seq)):
        if i in plot:
            annotate += '*'
        else:
            annotate += '-'
    print annotate
    print '\n'
    if '*' in annotate:
        return 1
    else:
        return 0

def find_tmembrane(records):
    """ identify transmembrane domains based on the following rules:
            (1) domain must either have 16 consecutive neutrally-charged residues or
            (2) if (1) is not met, allow for lysine (K) residues at locations n to n+2 and/or n+13 to n+15
    """

    num_caught = 0
    for rec in records:
        plot = []
        seq = records[rec].seq
        for i in range(3, len(seq)-12):
            if hydrophobicity(seq[i-3:i], seq[i:i+10], seq[i+10:i+13]):
                plot += range(i-3, i+13)

        print rec
        num_caught += print_seq(plot, seq)
    print num_caught
    print len(records)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='find phage transmembrane domains')
    parser.add_argument('gff3', type=file, help='GFF3 output of TMHMM')
    parser.add_argument('fasta', type=file, help='fasta file of protein(s)')
    args = parser.parse_args()
    taper_list(**vars(args))
