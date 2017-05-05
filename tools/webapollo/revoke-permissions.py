#!/usr/bin/env python
import time
import argparse
from webapollo import WebApolloInstance, AssertAdmin, AssertUser


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Admin Username')
    parser.add_argument('password', help='WA Admin Password')
    parser.add_argument('email', help='User Email')

    parser.add_argument('idList', type=argparse.FileType('r'), help='List of User IDs')
    parser.add_argument('--dry_run', action='store_true')

    args = parser.parse_args()
    wa = WebApolloInstance(args.apollo, args.username, args.password)

    gx_user = AssertAdmin(AssertUser(wa.users.loadUsers(email=args.email)))

    s464 = [x.strip() for x in args.idList.readlines()]

    print('# User ID\tEmail\tAction\tOrganism ID\tOrganism Name')
    users = wa.users.loadUsers()

    orgs = wa.organisms.findAllOrganisms()
    for user in users:
        if str(user.userId) in s464:
            for org in user.orgPerms():
                print('%s\t%s\t%s\t%s\t%s' % (user.userId, user.username, 'revoke', org['id'], org['organism']))
                if not args.dry_run:
                    wa.users.updateOrganismPermission(user, org['organism'],
                                                      administrate=False,
                                                      write=False, read=False,
                                                      export=False)
                    time.sleep(.5)
