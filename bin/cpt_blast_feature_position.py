#!/usr/bin/env python
import argparse
from Bio import SeqIO
import logging
logging.basicConfig(level=logging.INFO)

def get_id(feature=None, parent_prefix=None):
    result = ""
    if parent_prefix is not None:
        result += parent_prefix + '|'
    if 'locus_tag' in feature.qualifiers:
        result += feature.qualifiers['locus_tag'][0]
    elif 'gene' in feature.qualifiers:
        result += feature.qualifiers['gene'][0]
    elif 'product' in feature.qualifiers:
        result += feature.qualifiers['product'][0]
    else:
        result += '%s_%s_%s' % (feature.location.start, feature.location.end,
                                feature.location.strand)
    return result

def locate_hits(database, blast_results):
    hit_ids = [x.split()[1] for x in blast_results.readlines()]
    for record in SeqIO.parse(database, 'genbank'):
        for feature in record.features:
            if get_id(feature) in hit_ids:
                print '\t'.join([record.id, get_id(feature), feature.type,
                                 feature.location.start, feature.location.end])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Match blast hits to positions in Genbank file', epilog="")
    parser.add_argument('database', type=file, help='Genbank Database file')
    parser.add_argument('blast_results', type=file, help='12+ Column Blast Results')

    args = vars(parser.parse_args())
    locate_hits(**args)
