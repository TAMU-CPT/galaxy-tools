#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO)


def filter_blast(blast_results, top_n=5):
    id_counts = {}
    for line in blast_results:
        # evalue = line.split('\t')[10]
        id = line.split('\t')[0]
        if id in id_counts:
            id_counts[id] += 1
        else:
            id_counts[id] = 1

        if id_counts[id] <= top_n:
            print line


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter blast results')
    parser.add_argument('blast_results', type=file, help='Tabular Blast Results')
    parser.add_argument('top_n', type=int, help='Top N hits')

    args = parser.parse_args()
    filter_blast(**vars(args))
