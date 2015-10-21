#!/usr/bin/env python
# import sys
import re
import argparse
import hashlib
from BCBio import GFF
from Bio.Blast import NCBIXML
from gff3 import feature_lambda
from collections import OrderedDict

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
                    'sbject_range' : (hsp.sbjct_start, hsp.sbjct_end),
                    'query_range' : (hsp.query_start, hsp.query_end),
                    'name' : blast_record.query,
                    'iter_num' : iter_num
                })
        blast.append(blast_gene)
    return blast

def filter_lone_clusters(clusters):
    """ Removes all clusters with only one member and those with no hits """
    filtered_clusters = {}
    for key in clusters:
        if len(clusters[key]) > 1 and len(key) > 0:
            filtered_clusters[key] = clusters[key]
    return filtered_clusters

def test_true(feature, **kwargs):
    return True

def parse_gff(gff3):
    """ Extracts strand and start location to be used in cluster filtering """
    gff_info = {}
    for rec in GFF.parse(gff3):
        for feat in feature_lambda(
            rec.features,
            test_true,
            {},
            subfeatures = False
        ):
            gff_info[feat.id] = {'strand': feat.strand,'start': feat.location.start}

    gff_info = OrderedDict(sorted(gff_info.items(), key=lambda k: k[1]['start']))
    for i, feat_id in enumerate(gff_info):
        gff_info[feat_id].update({'index': i})

    return dict(gff_info)

def all_same(genes_list):
    """ Returns True if all gene names in cluster are identical """
    return all(gene['name'] == genes_list[0]['name'] for gene in genes_list)

def remove_duplicates(clusters):
    """ Removes clusters with multiple members but only one gene name """
    filtered_clusters = {}
    for key in clusters:
        if all_same(clusters[key]):
            continue
        else:
            filtered_clusters[key] = clusters[key]
    return filtered_clusters

class IntronFinder(object):
    """ IntronFinder objects are lists that contain a list of hits for every gene """
    def __init__(self, gff3, blastp):
        self.blast = []
        self.clusters = {}
        self.gff_info = {}

        self.gff_info = parse_gff(gff3)
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
        self.clusters = filter_lone_clusters(clusters)

    def check_strand(self):
        """ filters clusters for genes on the same strand """
        filtered_clusters = {}
        for key in self.clusters:
            pos_strand = []
            neg_strand = []
            for gene in self.clusters[key]:
                if self.gff_info[gene['name']]['strand'] == 1:
                    pos_strand.append(gene)
                else:
                    neg_strand.append(gene)
            if len(pos_strand) == 0 or len(neg_strand) == 0:
                filtered_clusters[key] = self.clusters[key]
            else:
                if len(pos_strand) > 1:
                    filtered_clusters[key+'_+1'] = pos_strand
                if len(neg_strand) > 1:
                    filtered_clusters[key+'_-1'] = neg_strand
        return filtered_clusters

    def check_gene_gap(self):
        filtered_clusters = {}
        for key in self.clusters:
            hits_lists = []
            gene_added = False
            for gene in self.clusters[key]:
                for hits in hits_lists:
                    for hit in hits:
                        if abs(self.gff_info[gene['name']]['index'] - self.gff_info[hit['name']]['index']) <= 10:
                            hits.append(gene)
                            gene_added = True
                            break
                if gene_added == False:
                    hits_lists.append([gene])

            for i, hits in enumerate(hits_lists):
                if len(hits) >= 2:
                    filtered_clusters[key+'_'+str(i)] = hits
        return remove_duplicates(filtered_clusters) # call remove_duplicates somewhere else?

    # maybe figure out how to merge with check_gene_gap?
    # def check_seq_gap():

    # could simply see if all sbject_ranges don't intersect, but can't. (right?)
    # def check_seq_overlap():

    # def modify_gff3(?):
    # """ merge 2 or more seq records into one """

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Intron detection')
    parser.add_argument('gff3', type=file, help='GFF3 gene calls')
    parser.add_argument('blastp', type=file, help='blast XML protein results')
    args = parser.parse_args()

    # create new IntronFinder object based on user input
    ifinder = IntronFinder(**vars(args))
    ifinder.create_clusters()
    ifinder.clusters = ifinder.check_strand()
    ifinder.clusters = ifinder.check_gene_gap()

    with open('out.txt', 'w') as handle:
        import pprint; pprint.pprint(ifinder.clusters)
