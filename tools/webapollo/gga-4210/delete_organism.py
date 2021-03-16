#!/usr/bin/env python
from __future__ import print_function

import argparse
import logging
import os

from apollo import accessible_organisms
from apollo.util import GuessOrg, OrgOrGuess

from arrow.apollo import get_apollo_instance

from webapollo import UserObj, handle_credentials

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def IsRemote():
    return 'GALAXY_SHARED_DIR' not in os.environ or len(os.environ['GALAXY_SHARED_DIR'].lower().strip()) == 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to completely delete an organism')
    parser.add_argument('email', help='User Email')
    OrgOrGuess(parser)

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
    org = wa.organisms.show_organism(org_cn)

    wa.organisms.delete_features(org['id'])

    if IsRemote():
        print(wa.remote.delete_organism(org['commonName']))
    else:
        print(wa.organisms.delete_organism(org['id']))
