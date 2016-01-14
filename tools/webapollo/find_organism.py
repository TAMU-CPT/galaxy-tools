#!/usr/bin/env python
import json
import argparse
from webapollo import WebApolloInstance

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Admin Username')
    parser.add_argument('password', help='WA Admin Password')

    parser.add_argument('--commonName', help='Common Name')

    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    data = wa.organisms.findAllOrganisms()
    if args.commonName is not None:
        data = [o for o in data if o['commonName'] == args.commonName]

    print json.dumps(data, indent=2)
