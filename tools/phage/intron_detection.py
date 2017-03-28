#!/usr/bin/env python
import sys
import re
import itertools
import argparse
import hashlib
import svgwrite
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
        for alignment in blast_record.alignments:
            hit_gis = alignment.hit_id + alignment.hit_def
            gi_nos = [str(gi) for gi in re.findall('(?<=gi\|)\d{9}', hit_gis)]
            for hsp in alignment.hsps:
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
                    'identity': hsp.identities,
                    'iter_num': iter_num
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
                gff_info[feat.id] = {
                    'strand': feat.strand,
                    'start': feat.location.start,
                    'loc': feat.location,
                    'feat': feat,
                }

    gff_info = OrderedDict(sorted(gff_info.items(), key=lambda k: k[1]['start']))
    for i, feat_id in enumerate(gff_info):
        gff_info[feat_id].update({'index': i})

    return dict(gff_info), _rec


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

        (self.gff_info, self.rec) = parse_gff(gff3)
        self.blast = parse_xml(blastp)

    def create_clusters(self):
        """ Finds 2 or more genes with matching hits """
        clusters = {}
        for gene in self.blast:
            for hit in gene:
                if ' ' in hit:
                    hit = hit[0:hit.index(' ')]

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
                    filtered_clusters[key + '_+1'] = pos_strand
                if len(neg_strand) > 1:
                    filtered_clusters[key + '_-1'] = neg_strand
        return filtered_clusters

    def check_gene_gap(self):
        filtered_clusters = {}
        for key in self.clusters:
            hits_lists = []
            gene_added = False
            for gene in self.clusters[key]:
                for hits in hits_lists:
                    for hit in hits:
                        if abs(self.gff_info[gene['name']]['index'] - self.gff_info[hit['name']]['index']) == 1:
                            hits.append(gene)
                            gene_added = True
                            break
                if not gene_added:
                    hits_lists.append([gene])

            for i, hits in enumerate(hits_lists):
                if len(hits) >= 2:
                    filtered_clusters[key + '_' + str(i)] = hits
        log.debug("check_gene_gap %s -> %s", len(self.clusters), len(filtered_clusters))
        return remove_duplicates(filtered_clusters)  # call remove_duplicates somewhere else?

    # maybe figure out how to merge with check_gene_gap?
    # def check_seq_gap():

    # also need a check for gap in sequence coverage?
    def check_seq_overlap(self):
        filtered_clusters = {}
        for key in self.clusters:
            add_cluster = True
            sbjct_ranges = []
            for gene in self.clusters[key]:
                sbjct_ranges.append(gene['sbjct_range'])

            combinations = list(itertools.combinations(sbjct_ranges, 2))

            for pair in combinations:
                if len(set(range(pair[0][0], pair[0][1])) &
                       set(range(pair[1][0], pair[1][1]))) > 0:
                    add_cluster = False
                    break
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

    def draw_genes(self, name):
        height = 200 * len(self.clusters)
        dwg = svgwrite.Drawing(filename=name, size=("1500px", "%spx" % height), debug=True)
        genes = dwg.add(dwg.g(id='genes', fill='white'))

        sbjct_y = 10
        query_x = 10
        for i, key in enumerate(self.clusters):
            log.info('Done with %s', i)
            for j, gene in enumerate(sorted(self.clusters[key],
                                            key=lambda k: self.gff_info[k['name']]['start'],
                                            reverse=True)):
                if j == 0:
                    genes.add(dwg.rect(insert=(10, sbjct_y), size=(gene['sbjct_length'], 20), fill='blue'))

                genes.add(dwg.rect(
                    insert=(query_x, sbjct_y + 80),
                    size=(gene['query_length'], 20),
                    fill='green'
                ))
                genes.add(dwg.text(gene['name'], insert=(query_x, sbjct_y + 95)))

                p1 = (gene['sbjct_range'][0] + 10, sbjct_y + 20)
                p2 = (gene['sbjct_range'][1] + 10, sbjct_y + 20)
                p3 = (gene['query_range'][1] + query_x, sbjct_y + 80)
                p4 = (gene['query_range'][0] + query_x, sbjct_y + 80)
                identity = float(gene['identity']) / gene['sbjct_length']
                genes.add(dwg.polyline([p1, p2, p3, p4], fill='red', opacity='%s' % identity))

                dwg.save()
                query_x += (gene['query_length'] + 10)

            query_x = 10

            if i == 0:
                sbjct_y += 190
            else:
                sbjct_y += 200

    def output_gff3(self, clusters):
        rec = copy.deepcopy(self.rec)
        rec.features = []
        for cluster_idx, cluster_id in enumerate(clusters):
            # Get the list of genes in this cluster
            associated_genes = set([x['name'] for x in clusters[cluster_id]])
            # Get the gene locations
            assoc_gene_info = {x: self.gff_info[x]['loc'] for x in associated_genes}
            # Now we construct a gene from the children as a "standard gene model" gene.

            # Get the minimum and maximum locations covered by all of the children genes
            gene_min = min([min(x[1].start, x[1].end) for x in assoc_gene_info.iteritems()])
            gene_max = max([max(x[1].start, x[1].end) for x in assoc_gene_info.iteritems()])

            evidence_notes = []
            for cluster_elem in clusters[cluster_id]:
                note = '{name} had {ident}% identity to GI:{pretty_gi}'.format(
                    pretty_gi=', '.join(cluster_elem['gi_nos']),
                    ident=int(100 * float(cluster_elem['identity']) / abs(cluster_elem['query_range'][1] - cluster_elem['query_range'][0])),
                    **cluster_elem
                )
                evidence_notes.append(note)

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
                    'ID': ['gp_%s.mRNA' % cluster_idx]
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
                # Set the qualifiers appropriately
                cds.qualifiers = {
                    'ID': ['gp_%s.CDS.%s' % (cluster_idx, idx)],
                    'score': score,
                }
                cdss.append(cds)

            # And we attach the things properly.
            mRNA.sub_features = cdss
            gene.sub_features = [mRNA]
            # And append to our record
            rec.features.append(gene)
        return rec


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Intron detection')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='GFF3 gene calls')
    parser.add_argument('blastp', type=argparse.FileType("r"), help='blast XML protein results')
    parser.add_argument('--svg', help='Path to output svg file to', default='clusters.svg')
    args = parser.parse_args()

    # create new IntronFinder object based on user input
    ifinder = IntronFinder(args.gff3, args.blastp)
    ifinder.create_clusters()
    ifinder.clusters = ifinder.check_strand()
    # ifinder.clusters = ifinder.check_gene_gap()
    # ifinder.clusters = ifinder.check_seq_overlap()

    condensed_report = ifinder.cluster_report()
    ifinder.draw_genes(args.svg)
    rec = ifinder.output_gff3(ifinder.clusters)
    GFF.write([rec], sys.stdout)

    # import pprint; pprint.pprint(ifinder.clusters)
