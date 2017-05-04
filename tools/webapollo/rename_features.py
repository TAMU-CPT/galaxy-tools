#!/usr/bin/env python
import argparse
import sys
from webapollo import WebApolloInstance
from webapollo import WAAuth, OrgOrGuess, GuessOrg, AssertUser, retry
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    WAAuth(parser)
    parser.add_argument('email', help='User Email')
    parser.add_argument('--set', choices=['name', 'qualifier'], help="Which attribute to set.")
    parser.add_argument('--qualifier_name', default="locus_tag")
    OrgOrGuess(parser)

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

    # Call setSequence to tell apollo which organism we're working with
    wa.annotations.setSequence(org['commonName'], org['id'])
    # Then get a list of features.
    features = wa.annotations.getFeatures()
    # For each feature in the features
    type_counts = {}
    sys.stdout.write("# Feature ID\tNew Name\n")
    for feature in sorted(features['features'], key=lambda x: x['location']['fmin']):
        if 'parent_type' in feature and feature['type']['name'] != 'tRNA':
            ft = feature['parent_type']['name']
        else:
            ft = feature['type']['name']

        if ft not in type_counts:
            type_counts[ft] = 0
        type_counts[ft] += 1

        new_name = '%s-%03d' % (ft, type_counts[ft])
        if args.set == 'name':
            def fn0():
                wa.annotations.setName(
                    feature['parent_id'],
                    new_name
                )
            def fn1():
                wa.annotations.setName(
                    feature['uniquename'],
                    new_name
                )
            retry(fn0)
            retry(fn1)
        else:
            def fn0():
                wa.annotations.addAttributes(
                    feature['parent_id'],
                    {'locus_tag': [new_name]}
                )
            retry(fn0)
        sys.stdout.write("%s\t%s\n" % (feature['uniquename'], new_name))
