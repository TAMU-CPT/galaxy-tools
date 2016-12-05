#!/usr/bin/env python
import argparse
import logging
from webapollo import WebApolloInstance
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Admin Username')
    parser.add_argument('password', help='WA Admin Password')

    parser.add_argument('cn', help='Organism Common Name')
    parser.add_argument('gaf_file', type=argparse.FileType("r"), help='GAF file')

    parser.add_argument('--email')

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

    # [(0, 'CPT'), (1, 'UUID'), (2, 'Unk'), (3, ''), (4, 'GO:0006260'), (5, 'GOA:interpro|GO_REF:0000002'), (6, 'IEA'), (7, 'InterPro:IPR008770'), (8, 'B'), (9, 'Phi-29 DNA terminal protein GP3'), (10, ''), (11, 'protein'), (12, 'taxon:0'), (13, '20160128'), (14, 'Pfam'), (15, ''), (16, '\n')]

    for line in args.gaf_file:
        if line.startswith('!'):
            continue
        data = line.split('\t')
        # Must be a CPT GAF annotation
        # TODO: eventually we'll want to push to CACAO
        assert data[0] == 'CPT'

        print list(enumerate(data))
