#!/usr/bin/env python
import logging
import copy
import argparse
import StringIO
import hashlib
from Bio import SeqIO
logging.basicConfig(level=logging.INFO)


def dedup(fasta_file, mutation='mutate'):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    output = StringIO.StringIO()

    known_records = {}
    ordering_keys = []

    for record in records:
        seq = str(record.seq).upper()
        md5 = hashlib.md5(seq).hexdigest()

        if md5 in known_records:
            known_records[md5]['keys'].append(record.id)
        else:
            known_records[md5] = {
                'keys': [],
                'rec': copy.deepcopy(record),
            }
            ordering_keys.append(md5)

    for key in ordering_keys:
        if len(known_records[key]['keys']) > 0:
            ident_str = ', '.join(known_records[key]['keys'])
            known_records[key]['rec'].description += " [Ident: %s]" % ident_str
        SeqIO.write(known_records[key]['rec'], output, "fasta")
    print output.getvalue()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='remove duplicate sequences in a fasta file')
    parser.add_argument('fasta_file', type=file, help='Fasta file')

    args = parser.parse_args()
    dedup(**vars(args))
