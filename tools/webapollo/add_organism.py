#!/usr/bin/env python
import json
import argparse
from webapollo import WebApolloInstance

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Admin Username')
    parser.add_argument('password', help='WA Admin Password')

    parser.add_argument('cn', help='Organism Common Name')
    parser.add_argument('jbrowse', help='JBrowse Data Directory')
    parser.add_argument('email', help='User Email')
    parser.add_argument('--blatdb', help='BlatDB Directory')
    parser.add_argument('--genus', help='Organism Genus')
    parser.add_argument('--species', help='Organism Species')
    parser.add_argument('--public', action='store_true', help='Make organism public')

    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    orgs = wa.organisms.addOrganism(
        args.cn,
        args.jbrowse,
        blatdb=args.blatdb,
        genus=args.genus,
        species=args.species,
        public=args.public
    )
    print json.dumps([org for org in orgs if org['commonName'] == args.cn], indent=2)
