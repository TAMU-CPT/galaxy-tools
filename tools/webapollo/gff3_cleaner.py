#!/usr/bin/env python
import sys
import argparse
from gff3 import feature_lambda, feature_test_type
from BCBio import GFF
import logging
logging.basicConfig(level=logging.WARN)
log = logging.getLogger(name='pav')

def coding_genes(feature_list):
    for x in feature_lambda(feature_list, feature_test_type, {'type': 'gene'}, subfeatures=True):
        if len(list(feature_lambda(x.sub_features, feature_test_type, {'type': 'CDS'}, subfeatures=False))) > 0:
            yield x


def genes(feature_list, feature_type='gene'):
    for x in feature_lambda(feature_list, feature_test_type,
                            {'type': feature_type},
                            subfeatures=True):
        yield x


def fix_apollo_issues(annotations, user_email):
    for rec in GFF.parse(annotations):
        for feat in rec.features:
            if feat.type != 'gene':
                continue

            for sf in feat.sub_features:
                if sf.type != 'mRNA':
                    continue

                for ssf in sf.sub_features:
                    if ssf.type != 'exon':
                        continue

                    if len(ssf) > 10:
                        continue

                    ssf.type = 'Shine_Dalgarno_sequence'

                sf.sub_features = [x for x in sf.sub_features if x.type not in
                                   ('non_canonical_five_prime_splice_site',
                                    'non_canonical_three_prime_splice_site')]
        yield rec

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('annotations', type=file, help='Parent GFF3 annotations')
    # parser.add_argument('genome', type=file, help='Genome Sequence')
    parser.add_argument('--user_email')

    args = parser.parse_args()
    for rec in fix_apollo_issues(**vars(args)):
        rec.annotations = {}
        GFF.write([rec], sys.stdout)
