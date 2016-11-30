#!/usr/bin/env python
import argparse
from BCBio import GFF
from Bio import SeqIO
import itertools


def taper_list(gff3, fasta):
    """ deletes records that already have identified transmembrane domains """

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
            for j in beg + end:
                if j not in aas and j != 'K':
                    return False
    return True


def ranges(i):
    """ makes a range out of a list of indices """

    for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
        b = list(b)
        yield b[0][1], b[-1][1]


def print_seq(locations, record):
    """ prints output """

    annotate = ""
    for i in range(len(record.seq)):
        if i in locations:
            annotate += '*'
        else:
            annotate += '-'

    if '*' in annotate:
        print record.id
        print record.seq
        print annotate
        for r in list(ranges(locations)):
            print r[0], '-', r[1]
        print '\n'
    else:
        return (record.id, record.seq)


def find_tmembrane(records):
    """ identify transmembrane domains based on the following rules:
            (1) domain must either have 16 consecutive neutrally-charged residues or
            (2) if (1) is not met, allow for lysine (K) residues at locations n to n+2 and/or n+13 to n+15
    """

    no_tmembrane_domains = []
    for rec in records:
        locations = []  # indices of hydrophobic domains
        for i in range(3, len(records[rec].seq) - 12):
            if hydrophobicity(records[rec].seq[i - 3:i], records[rec].seq[i:i + 10], records[rec].seq[i + 10:i + 13]):
                locations += [loc for loc in range(i - 3, i + 13) if loc not in locations]

        no_tmembrane_domains += [rec_id for rec_id in [print_seq(locations, records[rec])] if rec_id]

    print "Records with no found transmembrane domains:"
    for i in no_tmembrane_domains:
        print i[0]
        print i[1]
        print '\n'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='find phage transmembrane domains')
    parser.add_argument('gff3', type=file, help='GFF3 output of TMHMM')
    parser.add_argument('fasta', type=file, help='fasta file of protein(s)')
    args = parser.parse_args()
    taper_list(**vars(args))
