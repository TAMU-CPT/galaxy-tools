#!/usr/bin/env python
import logging
import random
import copy
logging.basicConfig(level=logging.INFO)
import argparse
from Bio import SeqIO, Seq
import StringIO


def mutate(char):
    targets = ['A', 'C', 'T', 'G']
    del targets[targets.index(char)]
    return targets

def insert(char):
    return [char + x for x in ('A', 'C', 'T', 'G')]

def delete(char):
    return [""]

def snp(fasta_file, mutation='mutate'):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    output = StringIO.StringIO()

    methodToCall = globals()[mutation]

    for record in records:
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
                    results.append(seq[0:i] + mutated + seq[i+1:])
                    messages.append('[%s at %s: %s]' % (mutation, i + 1, mutated))
                    ids.append('_%s%s%s' % (mutation[0], i + 1, mutated))

            if results is not None:
                for (result, message, id_add) in zip(results, messages, ids):
                    rec_copy = copy.deepcopy(record)
                    rec_copy.seq = Seq.Seq(result)
                    rec_copy.description += message
                    rec_copy.id += id_add
                    SeqIO.write(rec_copy, output, "fasta")
    print output.getvalue()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='make SNPs in fasta files')
    parser.add_argument('fasta_file', type=file, help='Fasta file')
    parser.add_argument('mutation', type=str, help='Type of mutation to make', default='mutate', choices=['mutate', 'insert', 'delete'])

    args = parser.parse_args()
    snp(**vars(args))
