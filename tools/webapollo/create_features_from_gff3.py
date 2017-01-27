#!/usr/bin/env python
import argparse
import json
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
    wa.annotations.setSequence(org['commonName'], org['id'])

    # print(wa.annotations.getFeatures())
    for rec in GFF.parse(args.gff3):
        for feature in rec.features:
            # featureData = featuresToFeatureSchema([feature])
            mRNA_data = featuresToFeatureSchema([feature.sub_features[0]])
            # print(json.dumps(mRNA_data, indent=2))
            print(wa.annotations.addTranscript(
                {
                    'features': mRNA_data
                }, trustme=True
            ))
