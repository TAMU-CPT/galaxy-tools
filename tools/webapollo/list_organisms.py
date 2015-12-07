#!/usr/bin/env python
import json
import argparse
from webapollo import WebApolloInstance

if __name__ == '__main__':
    json
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Username')
    parser.add_argument('password', help='WA Password')

    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    print json.dumps(wa.organisms.findAllOrganisms(), indent=2)

