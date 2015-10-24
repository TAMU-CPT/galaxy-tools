#!/usr/bin/env python
import sys
import os
from Bio import SeqIO


with open(sys.argv[1], 'r') as table_handle, open(sys.argv[2], 'r') as fasta_handle:
    data = []
    for line in table_handle:
        if line.startswith('#'):
            continue

        x = line.strip().split('\t')
        data.append({
            'round': x[0],
            'bin': x[1],
            'length': x[2],
            'fids': x[3].split(','),
        })

    def find_by_id(id):
        for x in data:
            if id in x['fids']:
                return x
        raise Exception("Unknown")

    os.makedirs('output')

    for record in SeqIO.parse(fasta_handle, 'fasta'):
        xd = find_by_id(record.id)
        fn = os.path.join('output', 'round_{round}_bin_{bin}.fasta'.format(**xd))
        with open(fn, 'a') as outhandle:
            SeqIO.write([record], outhandle, 'fasta')
