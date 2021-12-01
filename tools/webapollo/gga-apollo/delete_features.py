#!/usr/bin/env python
from __future__ import print_function

import argparse
import logging
import random
import json

from apollo import accessible_organisms
from apollo.util import GuessOrg, OrgOrGuess, retry

from arrow.apollo import get_apollo_instance

from webapollo import UserObj, handle_credentials

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to delete all features from an organism')
    parser.add_argument('email', help='User Email')
    parser.add_argument('--type', help='Feature type filter')
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
            args.org_id = x.strip()
            res = GuessOrg(args, wa)
            if res:
               org_cn.append(res[0])

    if len(org_cn) == 0:
        raise Exception("Organism Common Name not provided")


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

      sequences = wa.organisms.get_sequences(org['id'])
      for sequence in sequences['sequences']:
        log.info("Processing %s %s", org['commonName'], sequence['name'])
        # Call setSequence to tell apollo which organism we're working with
        wa.annotations.set_sequence(org['id'], sequence['name'])
        # Then get a list of features.
        features = wa.annotations.get_features()
        # For each feature in the features
        for feature in sorted(features['features'], key=lambda x: random.random()):
            if args.type:
                if args.type == 'tRNA':
                    if feature['type']['name'] != 'tRNA':
                        continue

                elif args.type == 'terminator':
                    if feature['type']['name'] != 'terminator':
                        continue

                elif args.type == 'mRNA':
                    if feature['type']['name'] != 'mRNA':
                        continue

                else:
                    raise Exception("Unknown type")

            # We see that deleteFeatures wants a uniqueName, and so we pass
            # is the uniquename field in the feature.
            def fn():
                wa.annotations.delete_feature(feature['uniquename'])
                print('Deleted %s [type=%s]' % (feature['uniquename'], feature['type']['name']))

            if not retry(fn, limit=3):
                print('Error %s' % feature['uniquename'])
