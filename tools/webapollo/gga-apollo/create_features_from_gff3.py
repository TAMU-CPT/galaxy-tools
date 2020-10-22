#!/usr/bin/env python
import argparse
import json
import logging

from apollo import accessible_organisms
from apollo.util import GuessOrg, OrgOrGuess

from arrow.apollo import get_apollo_instance

from webapollo import UserObj, handle_credentials
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('email', help='User Email')
    parser.add_argument('--source', help='URL where the input dataset can be found.')
    parser.add_argument('--use_name', action='store_true', help='Use the given name instead of generating one.')
    parser.add_argument('--disable_cds_recalculation', action='store_true', help='Disable CDS recalculation and instead use the one provided.')
    OrgOrGuess(parser)

    parser.add_argument('gff3', type=argparse.FileType('r'), help='GFF3 file')
    args = parser.parse_args()

    wa = get_apollo_instance()
    # User must have an account
    gx_user = UserObj(**wa.users._assert_or_create_user(args.email))
    handle_credentials(gx_user)

    # Get organism
    org_cn = GuessOrg(args, wa)
    if isinstance(org_cn, list):
        org_cn = org_cn[0]

    all_orgs = wa.organisms.get_organisms()
    if 'error' in all_orgs:
        all_orgs = []
    all_orgs = [org['commonName'] for org in all_orgs]
    if org_cn not in all_orgs:
        raise Exception("Could not find organism %s" % org_cn)

    orgs = accessible_organisms(gx_user, [org_cn], 'WRITE')
    if not orgs:
        raise Exception("You do not have write permission on this organism")

    load_result = wa.annotations.load_legacy_gff3(org_cn, args.gff3) # use_name=args.use_name, disable_cds_recalculation=args.disable_cds_recalculation)
    #print(json.dumps(load_result, indent=2))
