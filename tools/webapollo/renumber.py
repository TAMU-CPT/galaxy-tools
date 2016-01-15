#!/usr/bin/env python
import json
import copy
import argparse
from webapollo import WebApolloInstance
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Admin Username')
    parser.add_argument('password', help='WA Admin Password')

    parser.add_argument('cn', help='Organism Common Name')

    parser.add_argument('--email')
    parser.add_argument('--prefix', default='gene_')
    parser.add_argument('--leading', default='3')

    args = parser.parse_args()


    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    gx_user = wa.users.loadUsers(email=args.email)
    if len(gx_user) == 0:
        raise Exception("Unknown user. Please register first")

    # Must find the organism
    org = wa.organisms.findOrganismByCn(args.cn)
    # TODO: verify user has permissions on the organism
    wa.annotations.setSequence(args.cn, org['id'])
    raw_data =  wa.annotations.getFeatures()['features']

    data = sorted([
        (x['parent_id'], x['uniquename'], x['location']['fmin'], x['name'])
        for x in raw_data
    ], key=lambda x: x[2])

    format_string = args.prefix + '%0' + args.leading + 'd'
    format_string_mrna = format_string + '.mRNA'

    outData = copy.copy(org)
    outData['changes'] = []

    for i, feat in enumerate(data):
        idx = i + 1
        log.info('Renaming %s to %s', feat[3], format_string % idx)
        outData['changes'].append((feat[0], format_string % idx))
        wa.annotations.setName(feat[0], format_string % idx)
        outData['changes'].append((feat[1], format_string_mrna % idx))
        wa.annotations.setName(feat[1], format_string_mrna % idx)

    print json.dumps(outData, indent=2)
