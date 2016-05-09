#!/usr/bin/env python
# vim: set fileencoding=utf-8
import os
import json
import math
import argparse
import itertools
from gff3 import feature_lambda, feature_test_type, feature_test_quals, \
    coding_genes, genes, get_gff3_id, feature_test_location, get_rbs_from
from shinefind import NaiveSDCaller
from BCBio import GFF
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import reverse_complement, translate
from Bio.SeqFeature import SeqFeature, FeatureLocation
from jinja2 import Environment, FileSystemLoader
import itertools
from cpt import OrfFinder
import re
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name='pav')

# Path to script, required because of Galaxy.
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
# Path to the HTML template for the report

def annotation_table_report(record, wanted_cols):
    if wanted_cols is None or len(wanted_cols.strip()) == 0:
        return [], []

    def id(record, feature):
        """ID
        """
        return feature.id

    def name(record, feature):
        """Name
        """
        return feature.qualifiers.get('Name', ['None'])[0]

    def location(record, feature):
        """Location
        """
        return '{0.start}..{0.end}'.format(feature.location)

    def length(record, feature):
        """Length (AA)
        """
        cdss = list(genes(feature.sub_features, feature_type='CDS', sort=True))
        return str(sum([len(cds) for cds in cdss]) / 3)

    def notes(record, feature):
        """User entered Notes"""
        return feature.qualifiers.get('Note', [])

    def date_created(record, feature):
        """Created"""
        return feature.qualifiers.get('date_creation', ['None'])[0]

    def date_last_modified(record, feature):
        """Last Modified"""
        return feature.qualifiers.get('date_last_modified', ['None'])[0]

    def description(record, feature):
        """Description"""
        return feature.qualifiers.get('description', ['None'])[0]

    def owner(record, feature):
        """Owner

        User who created the feature. In a 464 scenario this may be one of
        the TAs."""
        return feature.qualifiers.get('owner', ['None'])[0]

    def product(record, feature):
        """Product

        User entered product qualifier (collects "Product" and "product"
        entries)"""
        return feature.qualifiers.get('product', feature.qualifiers.get('Product', ['None']))[0]

    def strand(record, feature):
        """Strand
        """
        return '+' if feature.location.strand > 0 else '-'

    def sd_spacing(record, feature):
        """Shine-Dalgarno spacing
        """
        rbss = get_rbs_from(gene)
        if len(rbss) == 0:
            return 'None'
        else:
            resp = []
            for rbs in rbss:
                cdss = list(genes(feature.sub_features, feature_type='CDS', sort=True))

                if rbs.location.strand > 0:
                    distance = min(cdss, key=lambda x: x.location.start - rbs.location.end)
                    distance_val = str(distance.location.start - rbs.location.end)
                    resp.append(distance_val)
                else:
                    distance = min(cdss, key=lambda x: x.location.end - rbs.location.start)
                    distance_val = str(rbs.location.start - distance.location.end)
                    resp.append(distance_val)

            if len(resp) == 1:
                return str(resp[0])
            return resp

    def sd_seq(record, feature):
        """Shine-Dalgarno sequence
        """
        rbss = get_rbs_from(gene)
        if len(rbss) == 0:
            return 'None'
        else:
            resp = []
            for rbs in rbss:
                resp.append(rbs.extract(record).seq)
            if len(resp) == 1:
                return str(resp[0])
            else:
                return resp

    def start_codon(record, feature):
        """Start Codon
        """
        cdss = list(genes(feature.sub_features, feature_type='CDS', sort=True))
        data = [x for x in cdss]
        if len(data) == 1:
            return str(data[0].extract(record).seq[0:3])
        else:
            return [
                '{0} ({1.location.start}..{1.location.end}:{1.location.strand})'.format(
                    x.extract(record).seq[0:3], x
                )
                for x in data
            ]

    def stop_codon(record, feature):
        """Stop Codon
        """
        return str(feature.extract(record).seq[-3:])

    def dbxrefs(record, feature):
        """DBxrefs
        """
        return feature.qualifiers.get('Dbxref', [])

    sorted_features = list(genes(record.features, sort=True))
    def upstream_feature(record, feature):
        """Next feature upstream"""
        upstream = None
        if feature.strand > 0:
            upstream_features = [x for x in sorted_features
                    if x.location.start < feature.location.start]
            if len(upstream_features) > 0:
                return upstream_features[-1]
            else:
                return None
        else:
            upstream_features = [x for x in sorted_features
                    if x.location.end > feature.location.end]

            if len(upstream_features) > 0:
                return upstream_features[0]
            else:
                return None

    def up_feat(record, feature):
        """Next feature upstream"""
        up = upstream_feature(record, feature)
        if up:
            return str(up)
        return 'None'


    def ig_dist(record, feature):
        """Distance to next feature on same strand"""
        up = upstream_feature(record, feature)
        if up:
            dist = None
            if feature.strand > 0:
                dist = feature.location.start - up.location.end
            else:
                dist = up.location.start - feature.location.end
            return str(dist)
        else:
            return 'None'


    cols = []
    data = []
    funcs = []
    lcl = locals()
    for x in [y.strip().lower() for y in wanted_cols.split(',')]:
        if x in lcl:
            funcs.append(lcl[x])
            # Keep track of docs
            func_doc = lcl[x].__doc__.strip().split('\n\n')
            # If there's a double newline, assume following text is the
            # "help" and the first part is the "name". Generate empty help
            # if not provided
            if len(func_doc) == 1:
                func_doc += ['']
            cols.append(func_doc)
        elif '__' in x:
           chosen_funcs = [lcl[y] for y in x.split('__')]
           func_doc = [' of '.join([y.__doc__.strip().split('\n\n')[0] for y in chosen_funcs[::-1]])]
           cols.append(func_doc)
           funcs.append(chosen_funcs)


    for gene in genes(record.features, sort=True):
        row = []
        for func in funcs:
            if isinstance(func, list):
                # If we have a list of functions, repeatedly apply them
                value = gene
                for f in func:
                    if value is None:
                        value = 'None'
                        break

                    value = f(record, value)
            else:
                # Otherwise just apply the lone function
                value = func(record, gene)

            if isinstance(value, list):
                value = [x.decode('utf-8') for x in value]
            else:
                value = value.decode('utf-8')

            row.append(value)
        # print row
        data.append(row)

    return data, cols

def evaluate_and_report(annotations, genome,
        reportTemplateName='phage_annotation_validator.html',
        annotationTableCols=''):
    """
    Generate our HTML evaluation of the genome
    """
    # Get features from GFF file
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    # Get the first GFF3 record
    # TODO: support multiple GFF3 files.
    record = list(GFF.parse(annotations, base_dict=seq_dict))[0]

    log.info("Producing an annotation table")
    annotation_table_data, annotation_table_col_names = annotation_table_report(record, annotationTableCols)

    # This is data that will go into our HTML template
    kwargs = {
        'annotation_table_data': annotation_table_data,
        'annotation_table_col_names': annotation_table_col_names,
    }

    def nice_strand(direction):
        if direction > 0:
            return '→'.decode('utf-8')
        else:
            return '←'.decode('utf-8')

    env = Environment(loader=FileSystemLoader(SCRIPT_PATH), trim_blocks=True, lstrip_blocks=True)
    env.filters['nice_id'] = get_gff3_id
    env.filters['nice_strand'] = nice_strand
    tpl = env.get_template(reportTemplateName)
    return tpl.render(**kwargs).encode('utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('annotations', type=file, help='Parent GFF3 annotations')
    parser.add_argument('genome', type=file, help='Genome Sequence')

    parser.add_argument('--reportTemplateName', help='Report template file name', default='phageqc_report_full.html')
    parser.add_argument('--annotationTableCols', help='Select columns to report in the annotation table output format')

    args = parser.parse_args()

    print evaluate_and_report(**vars(args))
