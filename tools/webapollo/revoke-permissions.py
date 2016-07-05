#!/usr/bin/env python
import time
import argparse
from webapollo import WebApolloInstance


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Admin Username')
    parser.add_argument('password', help='WA Admin Password')

    parser.add_argument('--commonName', help='Common Name')

    args = parser.parse_args()
    s464 = [
        42217,
        4637,
        4673,
        4793,
        4869,
        4904,
        4931,
        4958,
        4961,
        5036,
        5089,
        5092,
        5094,
        5116,
        5301,
        5427,
    ]

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    orgs = wa.organisms.findAllOrganisms()
    for org in orgs:
        for user in wa.users.loadUsers():
            if user.userId in s464:
                # interscetion of org + user
                print 'Revoking %s on %s' % (user, org['commonName'])
                wa.users.updateOrganismPermission(user, org['commonName'])
                time.sleep(.5)
