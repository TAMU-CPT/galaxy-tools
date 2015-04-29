#!/usr/bin/env python
import os
from Bio import SeqIO
import StringIO
import argparse
import logging
logging.basicConfig(level=logging.INFO)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta_file', type=file, help='Blastclust fasta input')
    parser.add_argument('blastclust', help='Blastclust Output')
    parser.add_argument('--ignore', type=int, help="Clusters smaller than this threshold will be ignored")
    args = parser.parse_args()

    id_to_clust = {}
    with open(args.blastclust, 'r') as handle:
        for i, line in enumerate(handle):
            for clust_id in line.strip().split():
                id_to_clust[clust_id] = i

    os.mkdir("out")
    for sequence in SeqIO.parse(args.fasta_file, 'fasta'):
        if sequence.id in id_to_clust:
            name = os.path.join('out', 'blastclust_%s.fa' % id_to_clust[sequence.id])
            with open(name, 'a') as handle:
                tmp = StringIO.StringIO("HI")
                SeqIO.write([sequence], tmp, 'fasta')
                handle.write(tmp.getvalue())
