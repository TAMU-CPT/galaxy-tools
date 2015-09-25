#!/usr/bin/env python
import sys
import re
import argparse
import hashlib
# from BCBio import GFF
from Bio.Blast import NCBIXML

def parse_xml(blastxml):
    """ Parses xml file to get desired info (genes, hits, etc) """
    blast = []
    for iter_num, blast_record in enumerate(NCBIXML.parse(blastxml)):
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
        blast.append(blast_gene)
    return blast

def filter_clusters(matches):
    """ Removes all clusters with only one member and those with no hits """
    filtered_matches = {}
    for key in matches:
        if len(matches[key]) > 1 and len(key) > 0:
            filtered_matches[key] = matches[key]
    return filtered_matches

class Hit(object):
    """ Each gene will have multiple Hits """
    def __init__(self, gi_nos, query_start, query_end, sbjct_start, sbjct_end, iter_def, iter_num):
        self.gi_nos = gi_nos
        self.query_start = query_start
        self.query_end = query_end
        self.sbjct_start = sbjct_start
        self.sbjct_end = sbjct_end
        self.iter_def = iter_def
        self.iter_num = iter_num

    def __repr__(self):
        return self.iter_def

class IntronFinder(object):
    """ IntronFinder objects are lists that contain a list of Hits for every gene """
    def __init__(self, genes, blastp):
        self.blast = []
        self.matches = {}
        self.blast = parse_xml(blastp)

    def create_clusters(self):
        """ Finds 2 or more genes with matching hits """
        clusters = {}
        for gene in self.blast:
            for hit in gene:
                name = hashlib.md5(','.join(hit.gi_nos)).hexdigest()
                if name in clusters:
                    if hit not in clusters[name]:
                        clusters[name].append(hit)
                else:
                    clusters[name] = [hit]
        self.matches = filter_clusters(clusters)

        # import pprint; pprint.pprint(self.matches)

    """ call checker functions to make sure intron is legit """
    # def check_strand():
    # def check_distance():

    # def modify_gff3(?):
    """ merge 2 or more seq records into one """

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Intron detection')
    parser.add_argument('genes', type=file, help='GFF3 Gene Calls')
    parser.add_argument('blastp', type=file, help='blast XML protein results')
    args = parser.parse_args()

    # create new IntronFinder object based on user input
    ifinder = IntronFinder(**vars(args))
    ifinder.create_clusters()

    with open('out.txt', 'w') as handle:
        import pprint; pprint.pprint(ifinder.matches)
