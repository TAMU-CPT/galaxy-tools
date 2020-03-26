#!/usr/bin/env python
import argparse
from webapollo import WebApolloInstance
from webapollo import WAAuth, OrgOrGuess, GuessOrgMulti, AssertUser, accessible_organisms
import logging

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
      assert x in orgs

    # The other person must already be an apollo user
    other_user = AssertUser(
        wa.users.loadUsers(email=args.share_with.replace("__at__", "@"))
    )

    # Did not appear referenced, suspect vestigial
    # org = wa.organisms.findOrganismByCn(org_cn)


    for x in org_cn:

      wa.users.updateOrganismPermission(
        other_user,
        x,
        administrate=False,
        write=args.write,
        export=args.export,
        read=args.read,
      )
