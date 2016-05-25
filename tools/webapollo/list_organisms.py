#!/usr/bin/env python
import json
import argparse
from webapollo import WAAuth, WebApolloInstance

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    print json.dumps(wa.organisms.findAllOrganisms(), indent=2)
