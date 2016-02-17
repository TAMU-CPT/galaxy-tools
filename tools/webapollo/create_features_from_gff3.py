#!/usr/bin/env python
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
    parser.add_argument('email', help='User Email')
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    gx_user = wa.users.loadUsers(email=args.email)
    if len(gx_user) == 0:
        raise Exception("Unknown user. Please register first")

    # TODO: Check user perms on org.

    org = wa.organisms.findOrganismByCn(args.cn)
    import pprint; pprint.pprint( org )

    wa.annotations.setSequence(args.cn, org['id'])
    from asdf import transcript, gene_transript
    wa.annotations.addFeature(
        {
            'features':[
                gene_transript
            ]
        }, trustme=True
    )

    # wa.annotations.addTranscript(
        # transcript, trustme=True
    # )
