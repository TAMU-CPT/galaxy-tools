#!/usr/bin/env python
from __future__ import print_function

import argparse
import logging
import sys
import json
from pathlib import Path

from apollo import accessible_organisms
from apollo.util import GuessOrg, OrgOrGuess

from arrow.apollo import get_apollo_instance

from webapollo import WebApolloInstance, WAAuth, AssertUser, UserObj, handle_credentials, galaxy_list_users

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sample script to delete all features from an organism"
    )
    WAAuth(parser)
    parser.add_argument("email", help="User Email")
    parser.add_argument("share_with", help="Share with this user (by email)")
    parser.add_argument("--write", action="store_true", help="Write permission")
    parser.add_argument("--export", action="store_true", help="Export permission")
    parser.add_argument("--read", action="store_true", help="Read permission")
    parser.add_argument('--org_json', type=argparse.FileType("r"), help='Apollo JSON output, source for common name')
    parser.add_argument('--org_raw', help='Common Name')
    parser.add_argument('--org_id', help='Organism ID')
    args = parser.parse_args()

    wa = get_apollo_instance()
    # User must have an account
    args.share_with = args.share_with.replace("__at__", "@")
    args.share_with = args.share_with.split(",")
    valUsers = []
    for users in args.share_with:
      other_user = wa.users.show_user(users.strip())
      if other_user == []:
        sys.stderr.write("Error: No such user " + users.strip() + " to share organism with, exiting...")
        exit(1)
      else:
        valUsers.append(other_user)

    orgs = wa.users.get_organism_permissions(args.email)
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

    # Get organism
    valOrgs = []

    for x in orgs:
        for y in org_cn:
          if y == str(x['organism']).strip() and 'EXPORT' in x['permissions']:
            valOrgs.append(y)

    for x in valOrgs:
      for users in valUsers:
        wa.users.update_organism_permissions(
          users['username'],
          x,
          administrate=False,
          write=args.write,
          export=args.export,
          read=args.read,
        )
        print("Successfully shared " + str(x) + " with user " + str(args.users['username']))
