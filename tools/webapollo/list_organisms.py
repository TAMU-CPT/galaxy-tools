#!/usr/bin/env python
import json
import argparse
from webapollo import WAAuth, WebApolloInstance, AssertUser

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)
    parser.add_argument('email', help='User Email')
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    gx_user = AssertUser(wa.users.loadUsers(email=args.email))
    # {u'organism': u'lso_phage_2', u'permissions': u'[]', u'userId': 142792},
    permissionMap = {
        x['organism']: x['permissions']
        for x in gx_user.organismPermissions
        if 'WRITE' in x['permissions'] or 'READ' in x['permissions']
    }

    orgs = wa.organisms.findAllOrganisms()
    orgs = [
        org for org in orgs
        if org['commonName'] in permissionMap
    ]

    print json.dumps(orgs, indent=2)
