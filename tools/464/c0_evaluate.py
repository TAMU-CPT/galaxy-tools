#!/usr/bin/env python
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
from guanine import GuanineClient
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
g = GuanineClient()


def validate(gff3, user_email):
    STUDENT_PASSES = False
    for rec in GFF.parse(gff3):
        for feature in feature_lambda(
            rec.features,
            feature_test_type,
            {'type': 'gene'},
            subfeatures=True,
        ):
            checks = []
            # Ownership
            if user_email not in feature.qualifiers.get('owner', []):
                continue
            else:
                checks.append(('ownership', True))

            checks.append((
                'color',
                '#500000' in feature.qualifiers.get('color', [])
            ))

            # dbxrefs
            checks.append((
                'pubmed',
                'PMID:26711672' in feature.qualifiers.get('Dbxref', [])
            ))
            checks.append((
                'go',
                'GO:0098009' in feature.qualifiers.get('Dbxref', [])
            ))
            checks.append((
                'dbxrefs',
                'CPT:283675' in feature.qualifiers.get('Dbxref', [])
            ))

            # Notes
            checks.append((
                'Note',
                'Howdy!' in feature.qualifiers.get('Note', [])
            ))

            # This gene passes
            if all([x[1] for x in checks]):
                STUDENT_PASSES = True

            for check in checks:
                print '%-10s %-10s %s' % (
                    feature.qualifiers.get('Name', ['unknown'])[0],
                    check[0],
                    'OK' if check[1] else 'NOT OK'
                )

    if STUDENT_PASSES:
        # Submit score to GUANINE
        g.submit(user_email, 'C0', 1.0)
    else:
        print "You have not successfully completed this exercise. Please try again. No annotations owned by you were found"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='verify against expectations')
    parser.add_argument('gff3', type=file, help='GFF3 annotations')
    parser.add_argument('user_email', help='User email')
    args = parser.parse_args()
    validate(**vars(args))
