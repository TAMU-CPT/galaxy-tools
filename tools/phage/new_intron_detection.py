#!/usr/bin/env python
import sys
import re
import argparse
# from BCBio import GFF
from Bio.Blast import NCBIXML

class Hit(object):

    def __init__(self, hits_list, query_start, query_end, sbjct_start, sbjct_end):
        self.hits_list = hits_list
        self.query_start = query_start
        self.query_end = query_end
        self.sbjct_start = sbjct_start
        self.sbjct_end = sbjct_end

    # def __repr__(self):
        # return '[%s], %d, %d, %d, %d' % (','.join(self.hits_list), self.query_start, self.query_end, self.sbjct_start, self.sbjct_end)

class IntronFinder(object):

    def __init__(self, genes, blastp):
        self.hits = []

        for blast_record in NCBIXML.parse(blastp):
            blast_iteration = []

            for alignment in blast_record.alignments:
                hit_gis = alignment.hit_id + alignment.hit_def
                hits_list =  [str(gi) for gi in re.findall('(?<=gi\|)\d{9}', hit_gis)]

                for hsp in alignment.hsps:
                    blast_iteration.append(Hit(hits_list,
                                               hsp.sbjct_start,
                                               hsp.sbjct_end,
                                               hsp.query_start,
                                               hsp.query_end))
            self.hits.append(blast_iteration)

    # checker function(s)
    # def check_strand():
    # def check_distance():

    # look through gi nos, find matching blast hits
    # call checker functions to make sure intron is legit
    # delete records from data structure that are not introns
    # def find_introns():

    # merge 2 or more seq records into one
    # def modify_gff3(?):

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Intron detection')
    parser.add_argument('genes', type=file, help='GFF3 Gene Calls')
    parser.add_argument('blastp', type=file, help='blast XML protein results')
    args = parser.parse_args()

    # create new IntronFinder objecti based on user input
    ifinder = IntronFinder(**vars(args))
