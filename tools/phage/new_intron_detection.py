#!/usr/bin/env python
# import sys
import re
import argparse
import hashlib
from BCBio import GFF
from Bio.Blast import NCBIXML
from gff3 import feature_lambda

def parse_xml(blastxml):
    """ Parses xml file to get desired info (genes, hits, etc) """
    blast = []
    for iter_num, blast_record in enumerate(NCBIXML.parse(blastxml), 1):
        blast_gene = []
        for alignment in blast_record.alignments:
            hit_gis = alignment.hit_id + alignment.hit_def
            gi_nos =  [str(gi) for gi in re.findall('(?<=gi\|)\d{9}', hit_gis)]
            for hsp in alignment.hsps:
                blast_gene.append({
                    'gi_nos' : gi_nos,
                    'sbject_start' : hsp.sbjct_start,
                    'sbjct_end' : hsp.sbjct_end,
                    'query_start' : hsp.query_start,
                    'query_end' : hsp.query_end,
                    'name' : blast_record.query,
                    'iter_num' : iter_num
                })
        blast.append(blast_gene)
    return blast

def filter_clusters(matches):
    """ Removes all clusters with only one member and those with no hits """
    filtered_matches = {}
    for key in matches:
        if len(matches[key]) > 1 and len(key) > 0:
            filtered_matches[key] = matches[key]
    return filtered_matches

def test_true(feature, **kwargs):
    return True

def strand(gff3):
    ids = {}
    for rec in GFF.parse(gff3):
        for feat in feature_lambda(
            rec.features,
            test_true,
            {},
            subfeatures = False
        ):
            ids[feat.id] = feat.strand
    return ids

class IntronFinder(object):
    """ IntronFinder objects are lists that contain a list of hits for every gene """
    def __init__(self, gff3, blastp):
        self.blast = []
        self.matches = {}
        self.lookup_strand = {}

        self.lookup_strand = strand(gff3)
        self.blast = parse_xml(blastp)

    def create_clusters(self):
        """ Finds 2 or more genes with matching hits """
        clusters = {}
        for gene in self.blast:
            for hit in gene:
                name = hashlib.md5(','.join(hit['gi_nos'])).hexdigest()
                if name in clusters:
                    if hit not in clusters[name]:
                        clusters[name].append(hit)
                else:
                    clusters[name] = [hit]
        self.matches = filter_clusters(clusters)


    def check_strand(self):
        """ filters clusters for genes on the same strand """
        refined_clusters = {}
        for key in self.matches:
            pos_strand = []
            neg_strand = []
            for gene in self.matches[key]:
                if self.lookup_strand[gene['name']] == 1:
                    pos_strand.append(gene)
                else:
                    neg_strand.append(gene)
            if len(pos_strand) == 0 or len(neg_strand) == 0:
                refined_clusters[key] = self.matches[key]
            else:
                if len(pos_strand) > 1:
                    refined_clusters[key+'_1'] = pos_strand
                if len(neg_strand) > 1:
                    refined_clusters[key+'_-1'] = neg_strand
        return refined_clusters

    # def check_distance():

    # def modify_gff3(?):
    """ merge 2 or more seq records into one """

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Intron detection')
    parser.add_argument('gff3', type=file, help='GFF3 gene calls')
    parser.add_argument('blastp', type=file, help='blast XML protein results')
    args = parser.parse_args()

    # create new IntronFinder object based on user input
    ifinder = IntronFinder(**vars(args))
    ifinder.create_clusters()
    ifinder.matches = ifinder.check_strand()

    with open('out.txt', 'w') as handle:
        import pprint; pprint.pprint(ifinder.matches)
