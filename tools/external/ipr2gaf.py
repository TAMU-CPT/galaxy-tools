#!/usr/bin/env python
import os
import argparse
import logging
from BCBio import GFF
from datetime import date
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
today = date.today()


def ipr2gaf(gff3):
    for record in GFF.parse(gff3):
        # Only want real features, not the fasta IPR sticks at the end
        if record.id.startswith('match$'):
            continue
        for feature in record.features:
            # Only want protein matches
            if feature.type != 'protein_match':
                continue

            base_columns = [
                # feature.qualifiers['source'][0],
                # feature.qualifiers['Name'][0],
                'CPT',
                record.id,
                'Unk',
                '',
            ]
            for go_term in feature.qualifiers.get('Ontology_term', []):
                for dbxref in feature.qualifiers['Dbxref']:
                    c = base_columns + [
                        go_term,
                        'GOA:interpro|GO_REF:0000002',
                        'IEA',
                        dbxref,
                        'B',
                        feature.qualifiers.get('signature_desc', [''])[0],
                        '',
                        'protein',
                        'taxon:0',
                        today.strftime("%Y%m%d"),
                        feature.qualifiers['source'][0],
                        '',
                        ''
                    ]
                    yield '\t'.join(c)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='IPR2GAF')
    parser.add_argument('gff3', help='GFF3 File')

    args = parser.parse_args()

    print '!gaf-version: 2.0'
    for gaf_line in ipr2gaf(args.gff3):
        print gaf_line
