#!/usr/bin/env python
import argparse
from webapollo import WebApolloInstance
from webapollo import WAAuth, AssertUser, retry


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='IRREVERSIBLY Delete organisms from Apollo')
    WAAuth(parser)
    parser.add_argument('email', help='User Email')
    parser.add_argument('orgids', help='Newline separated list of Organism IDs')

    args = parser.parse_args()

    # Process input file of organism ids
    with open(args.orgids, 'r') as f:
        orgids=[]
        for line in f:
            orgids.append(line.strip())

    # Sort and remove duplicate ids
    orgids = sorted(set(orgids))

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    gx_user = AssertUser(wa.users.loadUsers(email=args.email))

    # TODO: Check user perms on org.
    # finding all organisms first and the [x for x in all_organisms...] lines duplicate the action
    # of wa.organisms.FindOrganismByID() but reduces API calls, also allows us to skip incorrect IDs rather than break processing
    all_organisms = wa.organisms.findAllOrganisms()
    for orgid in orgids:
        orgs = [x for x in all_organisms if str(x['id']) == str(orgid)]
        if len(orgs) == 0:
            # Invalid ID number, skip it
            print("%s\t%s\t%s" % (orgid, "Not Found", "Not Found"))
            continue
        else:
            # ID was found, delete organism.
            org = orgs[0]
            #print("%s\t%s\t%s" % (org['id'], org['commonName'], "Deleted"))

            def fn():
                # Commenting the actual deletion for testing
                # wa.organisms.deleteOrganism(org['id'])
                print("%s\t%s\t%s" % (org['id'], org['commonName'], "Deleted"))

            if not retry(fn, limit=3):
                print("%s\t%s\t%s" % (org['id'], org['commonName'], "Error"))
