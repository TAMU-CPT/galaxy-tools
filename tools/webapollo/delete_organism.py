#!/usr/bin/env python
import argparse
from webapollo import WebApolloInstance
from webapollo import WAAuth, OrgOrGuess, GuessOrgMulti, AssertUser, accessible_organisms
import sys
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sample script to delete an organism"
    )
    WAAuth(parser)
    parser.add_argument("email", help="User Email")
    OrgOrGuess(parser)
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    gx_user = AssertUser(wa.users.loadUsers(email=args.email))
    

    # Get organism
    org_cn = GuessOrgMulti(args, wa)

    # Fetch all organisms
    all_orgs = wa.organisms.findAllOrganisms()
    # Figure out which are accessible to the user
    # filter out common names
    orgs = [d[0] for d in accessible_organisms(gx_user, all_orgs)]

    # This user MUST be allowed to access an organism before they can
    # modify permissions on it.
    for x in org_cn:
      sys.stderr.write("Assert " + str(x) + " is an accessible Apollo organism\n")
      assert x in orgs

    # "Remove" by setting all permissions to false, to preventusers who 
    # had the org shared with them from deleting it for everyone
    # Discuss if we want owner to be able to hard-delete, or if we want
    # to keep that ability in admin hands 
    for x in org_cn:
      wa.users.updateOrganismPermission(
        other_user,
        x,
        administrate=False,
        write=False,
        export=False,
        read=False,
      )
      print("Successfully removed " + str(org_cn) + " from organism list.")
    
