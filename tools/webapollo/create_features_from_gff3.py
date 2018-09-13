#!/usr/bin/env python
import sys
import time
import argparse
from webapollo import WebApolloInstance, featuresToFeatureSchema
from webapollo import WAAuth, OrgOrGuess, GuessOrg, AssertUser, retry
from BCBio import GFF
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)
    parser.add_argument('email', help='User Email')
    parser.add_argument('--source', help='URL where the input dataset can be found.')
    parser.add_argument('--name', default='product', help='Qualifer to use for the name of the added feature, default is \'product\'')
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

    bad_quals = ['date_creation', 'source', 'owner', 'date_last_modified', 'Name', 'ID']

    sys.stdout.write('# ')
    sys.stdout.write('\t'.join(['Feature ID', 'Apollo ID', 'Success', 'Messages']))
    sys.stdout.write('\n')
    # print(wa.annotations.getFeatures())
    for rec in GFF.parse(args.gff3):
        wa.annotations.setSequence(rec.id, org['id'])
        for feature in rec.features:
            # We can only handle genes right now
            if feature.type not in ('gene', 'terminator'):
                continue
            # Convert the feature into a presentation that Apollo will accept
            featureData = featuresToFeatureSchema([feature])
            if 'children' in featureData[0] and any([child['type']['name'] == 'tRNA' for child in featureData[0]['children']]):
                # We're experiencing a (transient?) problem where gene_001 to
                # gene_025 will be rejected. Thus, hardcode to a known working
                # gene name and update later.

                featureData[0]['name'] = 'tRNA_000'
                tRNA_sf = [child for child in feature.sub_features if child.type == 'tRNA'][0]
                tRNA_type = 'tRNA-' + tRNA_sf.qualifiers.get('Codon', ["Unk"])[0]

                if 'Name' in feature.qualifiers:
                    if feature.qualifiers['Name'][0].startswith('tRNA-'):
                        tRNA_type = feature.qualifiers['Name'][0]

                newfeature = wa.annotations.addFeature(featureData, trustme=True)

                def func0():
                    wa.annotations.setName(
                        newfeature['features'][0]['uniquename'],
                        tRNA_type,
                    )
                retry(func0)

                if args.source:
                    gene_id = newfeature['features'][0]['parent_id']

                    def setSource():
                        wa.annotations.addAttributes(gene_id, {'DatasetSource': [args.source]})

                    retry(setSource)

                sys.stdout.write('\t'.join([
                    feature.id,
                    newfeature['features'][0]['uniquename'],
                    'success',
                ]))
            elif featureData[0]['type']['name'] == 'terminator':
                # We're experiencing a (transient?) problem where gene_001 to
                # gene_025 will be rejected. Thus, hardcode to a known working
                # gene name and update later.
                featureData[0]['name'] = 'terminator_000'
                newfeature = wa.annotations.addFeature(featureData, trustme=True)
                newfeature = wa.annotations.addFeature(featureData, trustme=True)

                def func0():
                    wa.annotations.setName(
                        newfeature['features'][0]['uniquename'],
                        'terminator'
                    )

                retry(func0)

                if args.source:
                    newfeature2 = newfeature['features'][0]
                    gene_id = newfeature2.get('parent_id', newfeature2.get('uniquename', None))

                    def setSource():
                        wa.annotations.addAttributes(gene_id, {'DatasetSource': [args.source]})

                    retry(setSource)

                sys.stdout.write('\t'.join([
                    feature.id,
                    newfeature['features'][0]['uniquename'],
                    'success',
                ]))
            else:
                try:
                    # We're experiencing a (transient?) problem where gene_001 to
                    # gene_025 will be rejected. Thus, hardcode to a known working
                    # gene name and update later.
                    featureData[0]['name'] = 'gene_000'
                    # Extract CDS feature from the feature data, this will be used
                    # to set the CDS location correctly (apollo currently screwing
                    # this up (2.0.6))
                    CDS = featureData[0]['children'][0]['children']
                    CDS = [x for x in CDS if x['type']['name'] == 'CDS'][0]['location']
                    # Create the new feature
                    newfeature = wa.annotations.addFeature(featureData, trustme=True)
                    # Extract the UUIDs that apollo returns to us
                    mrna_id = newfeature['features'][0]['uniquename']
                    gene_id = newfeature['features'][0]['parent_id']
                    # Sleep to give it time to actually persist the feature. Apollo
                    # is terrible about writing + immediately reading back written
                    # data.
                    time.sleep(1)
                    # Correct the translation start, but with strand specific log
                    if CDS['strand'] == 1:
                        wa.annotations.setTranslationStart(mrna_id, min(CDS['fmin'], CDS['fmax']))
                    else:
                        wa.annotations.setTranslationStart(mrna_id, max(CDS['fmin'], CDS['fmax']) - 1)

                    # Finally we set the name, this should be correct.
                    time.sleep(0.5)
                    wa.annotations.setName(mrna_id, feature.qualifiers.get('product', feature.qualifiers.get('Name', ["Unknown"]))[0])
                    time.sleep(0.5)

                    def func():
                        wa.annotations.setName(gene_id, feature.qualifiers.get(args.name, feature.qualifiers.get('product', feature.qualifiers.get('Name', ["Unknown"])))[0])
                    retry(func)

                    if args.source:
                        gene_id = newfeature['features'][0]['parent_id']

                        def setSource():
                            wa.annotations.addAttributes(gene_id, {'DatasetSource': [args.source]})

                        retry(setSource)
                    extra_attr = {}
                    for (key, values) in feature.qualifiers.items():
                        if key in bad_quals:
                            continue

                        if key == 'Note':
                            def func2():
                                wa.annotations.addComments(gene_id, values)
                            retry(func2)
                        else:
                            extra_attr[key] = values

                    def func3():
                        wa.annotations.addAttributes(gene_id, extra_attr)
                    retry(func3)

                    sys.stdout.write('\t'.join([
                        feature.id,
                        gene_id,
                        'success',
                    ]))
                except Exception as e:
                    msg = str(e)
                    if '\n' in msg:
                        msg = msg[0:msg.index('\n')]
                    sys.stdout.write('\t'.join([
                        feature.id,
                        '',
                        'ERROR',
                        msg
                    ]))
            sys.stdout.write('\n')
