#!/usr/bin/env python
import json
import argparse
from webapollo import WAAuth, WebApolloInstance, AssertUser, accessible_organisms

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)
    parser.add_argument('email', help='User Email')
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    gx_user = AssertUser(wa.users.loadUsers(email=args.email))
    all_orgs = wa.organisms.findAllOrganisms()

    orgs = accessible_organisms(gx_user, all_orgs)

    print json.dumps(orgs, indent=2)
