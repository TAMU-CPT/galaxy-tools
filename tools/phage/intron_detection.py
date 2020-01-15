#!/usr/bin/env python
import sys
import re
import itertools
import argparse
import hashlib
import copy
from BCBio import GFF
from Bio.Blast import NCBIXML
from Bio.SeqFeature import SeqFeature, FeatureLocation
from gff3 import feature_lambda
from collections import OrderedDict
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()


def parse_xml(blastxml):
    """ Parses xml file to get desired info (genes, hits, etc) """
    blast = []
    discarded_records = 0
    for iter_num, blast_record in enumerate(NCBIXML.parse(blastxml), 1):
        blast_gene = []
        align_num = 0
        for alignment in blast_record.alignments:
            align_num += 1
            hit_gis = alignment.hit_id + alignment.hit_def
            gi_nos = [str(gi) for gi in re.findall('(?<=gi\|)\d{9,10}', hit_gis)]
            
            for hsp in alignment.hsps:
                #print(dir(hsp))
                x = float(hsp.identities) / (hsp.query_end - hsp.query_start)
                if x < .5:
                    discarded_records += 1
                    continue

                nice_name = blast_record.query 
                                             
                if ' ' in nice_name:
                    nice_name = nice_name[0:nice_name.index(' ')]

                blast_gene.append({
                    'gi_nos': gi_nos,
                    'sbjct_length': alignment.length,
                    'query_length': blast_record.query_length,
                    'sbjct_range': (hsp.sbjct_start, hsp.sbjct_end),
                    'query_range': (hsp.query_start, hsp.query_end),
                    'name': nice_name,
                    'evalue': hsp.expect,
                    'identity': hsp.identities,
                    'identity_percent': x,
                    'hit_num': align_num,
                    'iter_num': iter_num,
                    'match_id': alignment.title.partition(">")[0]
                })
        blast.append(blast_gene)
    log.debug("parse_blastxml %s -> %s", len(blast) + discarded_records, len(blast))
    return blast


def filter_lone_clusters(clusters):
    """ Removes all clusters with only one member and those with no hits """
    filtered_clusters = {}
    for key in clusters:
        if len(clusters[key]) > 1 and len(key) > 0:
            filtered_clusters[key] = clusters[key]
    log.debug("filter_lone_clusters %s -> %s", len(clusters), len(filtered_clusters))
    return filtered_clusters


def test_true(feature, **kwargs):
    return True


def parse_gff(gff3):
    """ Extracts strand and start location to be used in cluster filtering """
    log.debug("parse_gff3")
    gff_info = {}
    _rec = None
    for rec in GFF.parse(gff3):
        _rec = rec
        _rec.annotations = {}
        for feat in feature_lambda(
            rec.features,
            test_true,
            {},
            subfeatures=False
        ):
            if feat.type == 'CDS':
                if 'Name' in feat.qualifiers.keys():
                    CDSname = feat.qualifiers['Name']
                else:
                    CDSname = feat.qualifiers['ID']
                gff_info[feat.id] = {
                    'strand': feat.strand,
                    'start': feat.location.start,
                    'loc': feat.location,
                    'feat': feat,
                    'name': CDSname,
                }

   
    gff_info = OrderedDict(sorted(gff_info.items(), key=lambda k: k[1]['start']))
    endBase = 0
    for i, feat_id in enumerate(gff_info):
        gff_info[feat_id].update({'index': i})
        if gff_info[feat_id]['loc'].end > endBase:
            endBase = gff_info[feat_id]['loc'].end 
    

    return dict(gff_info), _rec, endBase


def all_same(genes_list):
    """ Returns True if all gene names in cluster are identical """
    return all(gene['name'] == genes_list[0]['name'] for gene in genes_list[1:])


def remove_duplicates(clusters):
    """ Removes clusters with multiple members but only one gene name """
    filtered_clusters = {}
    for key in clusters:
        if all_same(clusters[key]):
            continue
        else:
            filtered_clusters[key] = clusters[key]
    log.debug("remove_duplicates %s -> %s", len(clusters), len(filtered_clusters))
    return filtered_clusters


class IntronFinder(object):
    """ IntronFinder objects are lists that contain a list of hits for every gene """

    def __init__(self, gff3, blastp):
        self.blast = []
        self.clusters = {}
        self.gff_info = {}
        self.length = 0

        (self.gff_info, self.rec, self.length) = parse_gff(gff3)
        self.blast = parse_xml(blastp)

    def create_clusters(self):
        """ Finds 2 or more genes with matching hits """
        clusters = {}
        print(len(self.blast))
        for gene in self.blast:
            for hit in gene:
                if ' ' in hit['gi_nos']:
                    hit = hit[0:hit.index(' ')]
                name = hashlib.md5((','.join(hit['gi_nos'])).encode()).hexdigest()
                
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
                    filtered_clusters[key + '_+1'] = pos_strand
                if len(neg_strand) > 1:
                    filtered_clusters[key + '_-1'] = neg_strand
                    
        print(len(filtered_clusters))
        return filtered_clusters

    def check_gene_gap(self):
        filtered_clusters = {}
        for key in self.clusters:
            hits_lists = []
            gene_added = False
            for gene in self.clusters[key]:
                for hits in hits_lists:
                    for hit in hits:
                        if (abs(self.gff_info[gene['name']]['index'] - self.gff_info[hit['name']]['index']) <= 10) or ((len(self.gff_info) - (abs(self.gff_info[gene['name']]['index'] - self.gff_info[hit['name']]['index']))) <= 10): # Checks that they are within 10 array indices
                            hits.append(gene)  # of each other, including wrap around at the
                            gene_added = True  # end of the array.
                            break
                if not gene_added:
                    hits_lists.append([gene])

            for i, hits in enumerate(hits_lists):
                if len(hits) >= 2:
                    filtered_clusters[key + '_' + str(i)] = hits
        #for i in filtered_clusters:
         #   print(i)
          #  print(filtered_clusters[i])
        log.debug("check_gene_gap %s -> %s", len(self.clusters), len(filtered_clusters))
        return remove_duplicates(filtered_clusters)  # call remove_duplicates somewhere else?

    # maybe figure out how to merge with check_gene_gap?
    # def check_seq_gap():

    # also need a check for gap in sequence coverage?
    def check_seq_overlap(self, minimum = 0, maximum = 1000):
        filtered_clusters = {}
        for key in self.clusters:
            add_cluster = True
            sbjct_ranges = []
            for gene in self.clusters[key]:
                sbjct_ranges.append(gene['sbjct_range'])

            combinations = list(itertools.combinations(sbjct_ranges, 2))
            
            for pair in combinations:
                overlap = len(set(range(pair[0][0], pair[0][1])) & set(range(pair[1][0], pair[1][1])))
                minPair = pair[0]
                maxPair = pair[1]
                
                if minPair[0] > maxPair[0]:                
                  minPair = pair[1]
                  maxPair = pair[0]
                elif minPair[0] == maxPair[0] and minPair[1] > maxPair[1]:
                  minPair = pair[1]
                  maxPair = pair[0]
                if overlap > 0:
                  dist1 = maxPair[0] - minPair[0]
                  dist2 = max(maxPair[1], minPair[1]) - min(maxPair[1], minPair[1])
                else:  
                  dist1 = abs(maxPair[0] - minPair[1])
                  dist2 = (self.length - maxPair[1]) + minPair[0]  # Wraparound distance
                if minimum < 0:
                  if overlap > (minimum * -1):
                    add_cluster = False
                elif minimum == 0:
                  if overlap > 0:
                    add_cluster = False
                elif overlap > 0:
                  add_cluster = False

                if maximum < 0:
                  if overlap < (maximum * -1):
                    add_cluster = False
                elif maximum == 0:
                  if overlap == 0:
                    add_cluster = False

                if (dist1 > maximum or dist1 < minimum) and (dist2 > maximum or dist2 < minimum):
                  add_cluster = False
            if add_cluster:
                filtered_clusters[key] = self.clusters[key]
            
        log.debug("check_seq_overlap %s -> %s", len(self.clusters), len(filtered_clusters))
        return filtered_clusters

    def cluster_report(self):
        condensed_report = {}
        for key in self.clusters:
            for gene in self.clusters[key]:
                if gene['name'] in condensed_report:
                    condensed_report[gene['name']].append(gene['sbjct_range'])
                else:
                    condensed_report[gene['name']] = [gene['sbjct_range']]
        return condensed_report

    def cluster_report_2(self):
        condensed_report = {}
        for key in self.clusters:
            gene_names = []
            for gene in self.clusters[key]:
                gene_names.append((gene['name']).strip('CPT_phageK_'))
            if ', '.join(gene_names) in condensed_report:
                condensed_report[', '.join(gene_names)] += 1
            else:
                condensed_report[', '.join(gene_names)] = 1
        return condensed_report

    def cluster_report_3(self):
        condensed_report = {}
        for key in self.clusters:
            gene_names = []
            gi_nos = []
            for i, gene in enumerate(self.clusters[key]):
                if i == 0:
                    gi_nos = gene['gi_nos']
                gene_names.append((gene['name']).strip('.p01').strip('CPT_phageK_gp'))
            if ', '.join(gene_names) in condensed_report:
                condensed_report[', '.join(gene_names)].append(gi_nos)
            else:
                condensed_report[', '.join(gene_names)] = [gi_nos]
        return condensed_report

    
    def output_gff3(self, clusters):
        rec = copy.deepcopy(self.rec)
        rec.features = []
        for cluster_idx, cluster_id in enumerate(clusters):
            # Get the list of genes in this cluster
            associated_genes = set([x['name'] for x in clusters[cluster_id]])
            # print(associated_genes)
            # Get the gene locations
            assoc_gene_info = {x: self.gff_info[x]['loc'] for x in associated_genes}
            # Now we construct a gene from the children as a "standard gene model" gene.
            # Get the minimum and maximum locations covered by all of the children genes
            gene_min = min([min(x[1].start, x[1].end) for x in assoc_gene_info.items()])
            gene_max = max([max(x[1].start, x[1].end) for x in assoc_gene_info.items()])

            evidence_notes = []
            for cluster_elem in clusters[cluster_id]:
                note = '{name} had {ident}% identity to GI:{pretty_gi}'.format(
                    pretty_gi=', '.join(cluster_elem['gi_nos']),
                    ident=int(100 * float(cluster_elem['identity']) / abs(cluster_elem['query_range'][1] - cluster_elem['query_range'][0])),
                    **cluster_elem
                )
                evidence_notes.append(note)
            if gene_max - gene_min > .8 * float(self.length):
              evidence_notes.append("Intron is over 80% of the total length of the genome, possible wraparound scenario")
            # With that we can create the top level gene
            gene = SeqFeature(
                location=FeatureLocation(gene_min, gene_max),
                type='gene',
                id=cluster_id,
                qualifiers={
                    'ID': ['gp_%s' % cluster_idx],
                    'Notes': evidence_notes,
                }
            )
             
            # Below that we have an mRNA
            mRNA = SeqFeature(
                location=FeatureLocation(gene_min, gene_max),
                type='mRNA',
                id=cluster_id + '.mRNA',
                qualifiers={
                    'ID': ['gp_%s.mRNA' % cluster_idx],
                    'note': evidence_notes,
                }
            )

            # Now come the CDSs.
            cdss = []
            # We sort them just for kicks
            for idx, gene_name in enumerate(sorted(associated_genes, key=lambda x: int(self.gff_info[x]['start']))):
                # Copy the CDS so we don't muck up a good one
                cds = copy.copy(self.gff_info[gene_name]['feat'])
                # Get the associated cluster element (used in the Notes above)
                cluster_elem = [x for x in clusters[cluster_id] if x['name'] == gene_name][0]
                
                # Calculate %identity which we'll use to score
                score = int(1000 * float(cluster_elem['identity']) / abs(cluster_elem['query_range'][1] - cluster_elem['query_range'][0]))

                tempLoc = FeatureLocation(cds.location.start + (3 * (cluster_elem['query_range'][0])),
                                          cds.location.start + (3 * (cluster_elem['query_range'][1])),
                                          cds.location.strand)
                cds.location = tempLoc
                # Set the qualifiers appropriately
                cds.qualifiers = {
                    'ID': ['gp_%s.CDS.%s' % (cluster_idx, idx)],
                    'score': score,
                    'Name': self.gff_info[gene_name]['name'],
                    'evalue': cluster_elem['evalue'],
                    'Identity': cluster_elem['identity_percent'] * 100,
                    'intron_info': cluster_elem['match_id'], #'|'.join(cluster_elem['gi_nos']) + "| title goes here."
                }
                #cds.location.start = cds.location.start +
                cdss.append(cds)

            # And we attach the things properly.
            mRNA.sub_features = cdss
            gene.sub_features = [mRNA]
            # And append to our record
            rec.features.append(gene)
        return rec

    def output_xml(self, clusters):
        threeLevel = {}
        #print((clusters.viewkeys()))
        #print(type(enumerate(clusters)))
        #print(type(clusters))
        for cluster_idx, cluster_id in enumerate(clusters):
            #print(type(cluster_id))
            #print(type(cluster_idx)) 
            #print(type(clusters[cluster_id][0]['hit_num']))
            if not (clusters[cluster_id][0]['iter_num'] in threeLevel.keys):
                threeLevel[clusters[cluster_id][0]['iter_num']] = {}
        #for cluster_idx, cluster_id in enumerate(clusters):
        #    print(type(clusters[cluster_id]))
        #    b = {clusters[cluster_id][i]: clusters[cluster_id][i+1] for i in range(0, len(clusters[cluster_id]), 2)}
        #    print(type(b))#['name']))
        #for hspList in clusters:
        #for x, idx in (enumerate(clusters)):#for hsp in hspList:
        #    print("In X")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Intron detection')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='GFF3 gene calls')
    parser.add_argument('blastp', type=argparse.FileType("r"), help='blast XML protein results')
    parser.add_argument('--minimum', help='Gap minimum (Default 0, set to a negative number to allow overlap)', default = 0, type = int)
    parser.add_argument('--maximum', help='Gap minimum (Default 0, set to a negative number to allow overlap)', default = 1000, type = int)
    args = parser.parse_args()

    # create new IntronFinder object based on user input
    ifinder = IntronFinder(args.gff3, args.blastp)
    ifinder.create_clusters()
    ifinder.clusters = ifinder.check_strand()
    ifinder.clusters = ifinder.check_gene_gap()
    ifinder.clusters = ifinder.check_seq_overlap(minimum=args.minimum, maximum=args.maximum)
    #ifinder.output_xml(ifinder.clusters)
    #for x, idx in (enumerate(ifinder.clusters)):
    #print(ifinder.blast)

    condensed_report = ifinder.cluster_report()
    rec = ifinder.output_gff3(ifinder.clusters)
    GFF.write([rec], sys.stdout)
    
    

    #import pprint; pprint.pprint(ifinder.clusters)
