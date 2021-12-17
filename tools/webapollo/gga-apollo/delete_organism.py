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
    org_cn = []
    if args.org_json:
        org_cn = [x.get("commonName", None) for x in json.load(args.org_json)]
        org_cn = [x for x in orgs if x is not None]
    elif args.org_raw:
        args.org_raw = str(args.org_raw)
        args.org_raw = args.org_raw.split(",")
        for i in range(len(args.org_raw)):
            org_cn.append(args.org_raw[i].strip())
    elif args.org_id:
        orgList = str(args.org_id)
        orgList = orgList.split(",")
        for x in orgList:
            sys.stderr.write("Asserting " + x.strip() + " to be deleted")
            args.org_id = x.strip()
            res = GuessOrg(args, wa)
            if res:
               sys.stderr.write("Assert passed.")
               org_cn.append(res[0])

    if len(org_cn) == 0:
        raise Exception("Organism Common Name not provided")

    sys.stderr.write("CN List: " + str(org_cn))
            

    all_orgs = wa.organisms.get_organisms()
    if 'error' in all_orgs:
        all_orgs = []
    all_orgs = [org['commonName'] for org in all_orgs]
    
    for orgInd in org_cn:
      noOrg = True
      for eachOrg in all_orgs:
        if eachOrg == orgInd:
          noOrg = False
          break
      if noOrg:
        raise Exception("Could not find organism %s" % orgInd)
      else:
        orgs = accessible_organisms(gx_user, orgInd, 'WRITE')    
        if not orgs:
          raise Exception("You do not have write permission on this organism")

    for orgInd in org_cn:
      org = wa.organisms.show_organism(orgInd)

      wa.organisms.delete_features(org['id'])
 
      if IsRemote():
        print(wa.remote.delete_organism(org['commonName']))
      else:
        wa.organisms.delete_organism(org['id'], suppress_output=True)
        print(org)
