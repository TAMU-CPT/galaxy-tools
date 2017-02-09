#!/usr/bin/env python
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
# from guanine import GuanineClient


def validate(gff3):
    results = {}
    for rec in GFF.parse(gff3):
        for feature in feature_lambda(
            rec.features,
            feature_test_type,
            {'type': 'gene'},
            subfeatures=True,
        ):
            checks = []
            # dbxrefs
            checks.append((
                'dbxrefs',
                'CPT:283675' in feature.qualifiers.get('Dbxref', [])
            ))

            # Notes
            checks.append((
                'Note',
                'Howdy!' in feature.qualifiers.get('Note', [])
            ))

            owner = feature.qualifiers.get('owner', ['unknown'])[0]
            results[owner] = {
                'checks': checks,
                'score': 0,
            }

            # This gene passes
            results[owner]['score'] = [x[1] for x in checks].count(True)

    guanien_url = 'https://cpt.tamu.edu/guanine-backend/'
    token = auth('/galaxy/creds.json', guanine_url)
    for email, result in results.items():
        print(student, result)
        sid = student_id(email, guanine_url, token)
    # print(results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='verify against expectations')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='GFF3 annotations')
    args = parser.parse_args()
    validate(**vars(args))
