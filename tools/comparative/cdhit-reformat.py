#!/usr/bin/env python
import argparse
import logging
import os
from Bio import SeqIO
logging.basicConfig(level=logging.INFO)


def extract_clust(cdhit, min_size):
    cluster = []
    cluster_id = None
    output = 0
    skipped = 0

    for line in cdhit:
        if line.startswith('>'):
            if len(cluster) >= min_size:
                yield cluster, cluster_id
                output += 1
            else:
                if cluster_id:
                    skipped += 1
            cluster = []
            cluster_id = line.strip().split(' ')[1]
        elif cluster_id is not None:
            protein_id = line[line.index(', >') + 3:line.rindex('... ')]
            # Move representative to front.
            if '... *' in line:
                cluster = [protein_id] + cluster
            else:
                cluster.append(protein_id)

    if len(cluster) >= min_size:
        yield cluster, cluster_id
        output += 0
    else:
        skipped += 1
    print "Filtered %0.3f%% of clusters" % (100 * float(output) / (output + skipped))


def cdhit_process(cdhit, fasta, output_dir, min_size=1):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    for (cluster_elements, cluster_id) in extract_clust(cdhit, min_size):
        cluster_seqs = [seq_dict[x] for x in cluster_elements]
        with open(os.path.join(output_dir, cluster_id + '.fa'), 'w') as handle:
            SeqIO.write(cluster_seqs, handle, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cdhit', type=argparse.FileType("r"))
    parser.add_argument('fasta', type=argparse.FileType("r"))
    parser.add_argument('output_dir')
    parser.add_argument('--min_size', type=int, default=1)
    args = parser.parse_args()
    cdhit_process(**vars(args))
