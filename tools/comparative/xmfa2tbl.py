#!/usr/bin/env python
from Bio import SeqIO
from xmfa import parse_xmfa, percent_identity
import argparse
import logging
import itertools
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def _id_tn_dict(sequences):
    """Figure out sequence IDs AND sequence lengths from fasta file
    """
    label_convert = {}
    if sequences is not None:
        if len(sequences) == 1:
            for i, record in enumerate(SeqIO.parse(sequences[0], 'fasta')):
                label_convert[str(i + 1)] = {'id': record.id, 'len': len(record)}
        else:
            for i, sequence in enumerate(sequences):
                for record in SeqIO.parse(sequence, 'fasta'):
                    label_convert[str(i + 1)] = {'id': record.id, 'len': len(record)}
                    continue

    return label_convert


def total_similarity(xmfa_file, sequences=None, dice=False):
    if sequences is None:
        raise Exception("Must provide a non-zero number of sequence files")

    label_convert = _id_tn_dict(sequences)
    lcbs = parse_xmfa(xmfa_file)

    # make a matrix based on number of sequences
    table = [[0 for x in range(len(label_convert))] for y in range(len(label_convert))]

    for lcb in lcbs:
        # ignore LCBs containing only one sequence
        if len(lcb) == 0:
            continue

        # permutations based on num sequences to compare for current LCB
        compare_seqs = list(itertools.permutations(range(0, len(lcb)), 2))
        for permutation in compare_seqs:
            (i, j) = permutation
            similarity = percent_identity(lcb[i]['seq'], lcb[j]['seq'])
            # find length of sequence in LCB
            length_seq_lcb = lcb[i]['end'] - (lcb[i]['start'] - 1)
            # populate table with normalized similarity value based on length_seq_lcb
            table[int(lcb[i]['id']) - 1][int(lcb[j]['id']) - 1] += length_seq_lcb * similarity

    # finalize total percent similarity by dividing by length of parent sequence
    for i in range(len(label_convert)):
        for j in range(len(label_convert)):
            if dice:
                table[i][j] = 2 * table[i][j] / (label_convert[str(i + 1)]['len'] + label_convert[str(j + 1)]['len'])
            else:
                table[i][j] = table[i][j] / label_convert[str(i + 1)]['len']

    # insert 1 for comparisons between the same sequence
    for i in range(len(label_convert)):
        table[i][i] = 100

    # print table
    names = []
    for i in range(len(label_convert)):
        names.append(sorted(label_convert.values())[i]['id'])
    print '\t' + '\t'.join(names)
    for row in table:
        print names[table.index(row)] + '\t' + '\t'.join(map(str, row))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert XMFA alignments to gff3', prog='xmfa2gff3')
    parser.add_argument('xmfa_file', type=file, help='XMFA File')
    parser.add_argument('sequences', type=file, nargs='+',
                        help='Fasta files (in same order) passed to parent for reconstructing proper IDs')
    parser.add_argument('--dice', action='store_true', help='Use dice method for calculating % identity')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    total_similarity(**vars(args))
