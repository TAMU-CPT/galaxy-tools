#!/usr/bin/env python
# import sys
import re
import argparse
# from BCBio import GFF
from Bio.Blast import NCBIXML

class Hit(object):
    def __init__(self, gi_nos, query_start, query_end, sbjct_start, sbjct_end, iter_def, iter_num):
        self.gi_nos = gi_nos
        self.query_start = query_start
        self.query_end = query_end
        self.sbjct_start = sbjct_start
        self.sbjct_end = sbjct_end
        self.iter_def = iter_def
        self.iter_num = iter_num

    def __repr__(self):
        return '[%s], %d, %d, %d, %d' % (','.join(self.gi_nos), self.query_start, self.query_end, self.sbjct_start, self.sbjct_end)

class IntronFinder(object):
    def __init__(self, genes, blastp):
        self.blast = []
        self.matches = []

        for iter_num, blast_record in enumerate(NCBIXML.parse(blastp)):
            blast_gene = []
            for alignment in blast_record.alignments:
                hit_gis = alignment.hit_id + alignment.hit_def
                gi_nos =  [str(gi) for gi in re.findall('(?<=gi\|)\d{9}', hit_gis)]
                for hsp in alignment.hsps:
                    blast_gene.append(Hit(gi_nos,
                                          hsp.sbjct_start,
                                          hsp.sbjct_end,
                                          hsp.query_start,
                                          hsp.query_end,
                                          blast_record.query,
                                          iter_num))
            self.blast.append(blast_gene)

    # checker function(s)
    # call checker functions to make sure intron is legit
    # def check_strand():
    # def check_distance():

    # look through gi nos, find matching ones, append self.matches
    def compare_genes(self, gene1, gene2):
        for hit_1 in gene1:
            for hit_2 in gene2:
                # find intersection of hit_1.gi_nos and hit_2 gi_nos
                if len(set(hit_1.gi_nos) & set(hit_2.gi_nos)) > 0:
                    self.matches.append((hit_1, hit_2))

    # iterate through blast iterations, compare each gene to all other genes
    def find_introns(self):
        for i in range(0, len(self.blast)):
            print len(self.blast)
            print i
            for j in range(i, len(self.blast)):
                if i == j:
                    continue
                self.compare_genes(self.blast[i], self.blast[j])

       # test match list
        # with open('out.txt', 'w') as handle:
            # for tup in self.matches:
                # handle.write(tup[0].iter_def + ' ' + tup[1].iter_def)

    # merge 2 or more seq records into one
    # def modify_gff3(?):

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Intron detection')
    parser.add_argument('genes', type=file, help='GFF3 Gene Calls')
    parser.add_argument('blastp', type=file, help='blast XML protein results')
    args = parser.parse_args()

    # create new IntronFinder object based on user input
    ifinder = IntronFinder(**vars(args))
    ifinder.find_introns()
