#!/usr/bin/env python
from __future__ import print_function

import argparse
import logging
import os
import sys

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
    all_orgs = wa.organisms.get_organisms()
    if 'error' in all_orgs:
      raise Exception("Error while getting the list of organisms: %s" % all_orgs)
        

    org_cn = []
    #org_ids = []
    if args.org_json:
        org_cn = [x.get("commonName", None) for x in json.load(args.org_json)]
        org_cn = [x for x in orgs if x is not None]
    elif args.org_raw:
        args.org_raw = str(args.org_raw)
        args.org_raw = args.org_raw.split(",")
        for i in range(len(args.org_raw)):
            sys.stderr.write("Asserting " + args.org_raw[i].strip() + " to be deleted")
            org_cn.append(args.org_raw[i].strip())
    elif args.org_id:
        orgList = str(args.org_id)
        orgList = orgList.split(",")
        for orgInd in orgList:
          foundOrg = False
          orgInd = orgInd.strip()
          sys.stderr.write("Asserting " + orgInd.strip() + " to be deleted")
          for allOrgInd in all_orgs:
            if str(allOrgInd["id"]) == orgInd:
              org_cn.append(str(allOrgInd["commonName"]))
              foundOrg = True
              break        
          if not foundOrg:
            raise Exception("Organism '%s' not found" % orgInd)

    if len(org_cn) == 0:
        raise Exception("Organism Common Name not provided")
    
    sys.stderr.write("CN List: " + str(org_cn))
            

    # all_orgs = [org['commonName'] for org in all_orgs]
    permdOrgs = accessible_organisms(gx_user, org_cn, 'WRITE')    
    if len(permdOrgs) != len(org_cn):
          raise Exception("You do not have write permission on one or more of the supplied organisms")

    for orgInd in org_cn:
      org = wa.organisms.show_organism(orgInd)
      # Theoretically removable except direct entry would need a mechanism for getting ID in addition to common name

      wa.organisms.delete_features(org['id'])
 
      if IsRemote():
        print(wa.remote.delete_organism(org['commonName']))
      else:
        wa.organisms.delete_organism(org['id'], suppress_output=True)
        print(org)
