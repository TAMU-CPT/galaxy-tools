#!/usr/bin/env python
import argparse
import sys
import copy
from Bio import SeqIO
from Bio.Seq import Seq
import logging
logging.basicConfig(level=logging.INFO)


def mutate(char):
    targets = ['A', 'C', 'T', 'G']
    del targets[targets.index(char)]
    return targets


def insert(char):
    return [char + x for x in ('A', 'C', 'T', 'G')]


def delete(char):
    return [""]


def snp(fasta_file, mutation='mutate', translate=False):

    methodToCall = globals()[mutation]

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        # For every character
        for i in [-1] + range(len(seq)):
            results = []
            messages = []
            ids = []
            if i == -1:
                if mutation == 'insert':
                    for tmp in methodToCall(""):
                        results.append(tmp + seq)
                        messages.append('[%s at %s: %s]' % (mutation, i + 1, tmp))
                        ids.append('_%s%s%s' % (mutation[0], i + 1, tmp))
                else:
                    results = None
            else:
                # Get all possible mutations
                for mutated in methodToCall(seq[i]):
                    # And add those
                    results.append(seq[0:i] + mutated + seq[i + 1:])
                    messages.append('[%s at %s: %s]' % (mutation, i + 1, mutated))
                    ids.append('_%s%s%s' % (mutation[0], i + 1, mutated))

            if results is not None:
                seen = {}
                for (result, message, id_add) in zip(results, messages, ids):
                    rec_copy = copy.deepcopy(record)
                    if translate:
                        rec_copy.seq = Seq(result).translate(table=11)
                        y = str(rec_copy.seq)
                        if y in seen:
                            continue
                        else:
                            seen[y] = None
                    else:
                        rec_copy.seq = Seq(result)
                    rec_copy.description += message
                    rec_copy.id += id_add
                    yield [rec_copy]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate ALL possible SNPs in fasta sequences')
    parser.add_argument('fasta_file', type=argparse.FileType('r'), help='Fasta file')
    parser.add_argument('mutation', type=str, help='Type of mutation to make',
                        default='mutate', choices=['mutate', 'insert', 'delete'])
    parser.add_argument('--translate', action='store_true')

    args = parser.parse_args()
    for record in snp(**vars(args)):
        SeqIO.write(record, sys.stdout, "fasta")
