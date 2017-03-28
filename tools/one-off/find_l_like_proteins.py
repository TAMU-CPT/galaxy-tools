#!/usr/bin/env python
import argparse
from Bio import SeqIO


def find_hydrophobic_seq(seq):
    """
        Accomplishes #1 and #2 from find_l_like_proteins function
    """
    hydrophobic_residues = 'AGVLIPFHSTYCMWQN'  # allowed hydrophobic amino acids

    count = 0  # count for number of hydrophobic residues in a row
    start = 0
    end = 0
    domain = ''
    seq = seq + 'X'  # add a final non hydrophobic residue so loop will finish
    for num, s in enumerate(seq):
        if s in hydrophobic_residues:
            if not count:
                start = num
            count += 1
            end = num + 1
            domain = seq[start:end]
        else:
            if count >= 10:
                if domain.endswith('LS'):
                    yield {'start': start, 'end': end, 'domain': domain}
                elif 'LS' in domain and len(domain.rsplit('LS', 1)[0]) >= 8:
                    domain = domain.rsplit('LS', 1)[0] + 'LS'
                    end = start + len(domain)
                    yield {'start': start, 'end': end, 'domain': domain}
            count = 0


def n_terminus_charge(seq):
    """
        Accomplishes #3 from find_l_like_proteins function
    """

    positive_residues = 'RK'
    negative_residues = 'DE'

    charge = 0
    for letter in seq:
        if letter in positive_residues:
            charge += 1
        if letter in negative_residues:
            charge -= 1
    if charge >= 2:
        return True


def find_l_like_proteins(fasta):
    """ Returns proteins that have:
            1. >= 10 hydrophobic aa's in a row (AGVLIPFHSTYCMWQN) TODO: QN pref at the ends
            2. 'LS' right after that domain
            3. minimum net charge of +2 for N terminus
    """
    records = list(SeqIO.parse(fasta, "fasta"))

    for record in records:
        for a in find_hydrophobic_seq(record.seq):
            if n_terminus_charge(record.seq[0:a['start']]):
                if record.name.endswith('.CDS'):
                    print record.name[:-4]
                else:
                    print record.name


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Finds L-like proteins')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='Fasta file')

    args = parser.parse_args()
    find_l_like_proteins(args.fasta)
