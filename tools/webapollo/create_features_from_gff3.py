#!/usr/bin/env python
import sys
import json
import time
import argparse
from webapollo import WebApolloInstance, featuresToFeatureSchema
from webapollo import WAAuth, OrgOrGuess, GuessOrg, AssertUser
from BCBio import GFF
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)
    parser.add_argument('email', help='User Email')
    OrgOrGuess(parser)

    parser.add_argument('gff3', type=argparse.FileType('r'), help='GFF3 file')
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    gx_user = AssertUser(wa.users.loadUsers(email=args.email))

    # Get organism
    org_cn = GuessOrg(args, wa)
    if isinstance(org_cn, list):
        org_cn = org_cn[0]

    # TODO: Check user perms on org.
    org = wa.organisms.findOrganismByCn(org_cn)
    wa.annotations.setSequence('Mijalis', org['id'])

    bad_quals = ['date_creation', 'source', 'owner', 'date_last_modified', 'Name', 'ID']

    sys.stdout.write('# ')
    sys.stdout.write('\t'.join(['Feature ID', 'Apollo ID', 'Success', 'Messages']))
    sys.stdout.write('\n')

    # print(wa.annotations.getFeatures())
    for rec in GFF.parse(args.gff3):
        for feature in rec.features:
            featureData = featuresToFeatureSchema([feature])
            print featureData

            # try:
            featureData[0]['name'] = 'gene_000'
            CDS = featureData[0]['children'][0]['children']
            CDS = [x for x in CDS if x['type']['name'] == 'CDS'][0]['location']
            newfeature = wa.annotations.addFeature(featureData, trustme=True)
            mrna_id = newfeature['features'][0]['uniquename']
            gene_id = newfeature['features'][0]['parent_id']
            time.sleep(1)
            # Strand stuff
            if CDS['strand'] == 1:
                wa.annotations.setTranslationStart(mrna_id, min(CDS['fmin'], CDS['fmax']))
            else:
                wa.annotations.setTranslationStart(mrna_id, max(CDS['fmin'], CDS['fmax']) - 1)

            wa.annotations.setName(mrna_id, feature.qualifiers.get('product', ["Unknown"])[0])
            wa.annotations.setName(gene_id, feature.qualifiers.get('product', ["Unknown"])[0])
            # except Exception as e:
                # print e, feature.id
            # sys.exit()

            # try:
                # mRNA_data = featuresToFeatureSchema([feature.sub_features[0]])
                # resp = wa.annotations.addTranscript(
                    # {
                        # 'features': mRNA_data
                    # }, trustme=True
                # )

            # # Get the unique ID back from apollo
                # mRNA_uniq = resp['features'][0]['uniquename']
                # # Mapping table
                # sys.stdout.write('\t'.join([
                    # feature.id,
                    # mRNA_uniq,
                    # 'success',
                    # "Dropped qualifiers: %s" % (json.dumps({k: v for (k, v) in feature.qualifiers.items() if k not in bad_quals})),
                # ]))
            # except Exception as e:
                # sys.stdout.write('\t'.join([
                    # feature.id,
                    # '',
                    # 'UNKNOWN ERROR',
                    # str(e)
                # ]))

            # sys.stdout.write('\n')
