#!/usr/bin/env python
import sys
import json
import argparse
import time
from webapollo import WAAuth, WebApolloInstance, OrgOrGuess, GuessOrg, AssertUser
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create or update an organism in an Apollo instance')
    WAAuth(parser)

    parser.add_argument('jbrowse', help='JBrowse Data Directory')
    parser.add_argument('email', help='User Email')
    OrgOrGuess(parser)
    parser.add_argument('--genus', help='Organism Genus')
    parser.add_argument('--species', help='Organism Species')
    parser.add_argument('--public', action='store_true', help='Make organism public')
    # parser.add_argument('--group', help='Give access to a user group')

    args = parser.parse_args()
    wa = WebApolloInstance(args.apollo, args.username, args.password)

    org_cn = GuessOrg(args, wa)
    if isinstance(org_cn, list):
        org_cn = org_cn[0]

    # User must have an account
    gx_user = AssertUser(wa.users.loadUsers(email=args.email))

    log.info("Determining if add or update required")
    try:
        org = wa.organisms.findOrganismByCn(org_cn)
    except Exception:
        org = None

    if org:
        has_perms = False
        for user_owned_organism in gx_user.organismPermissions:
            if 'WRITE' in user_owned_organism['permissions']:
                has_perms = True
                break

        if not has_perms:
            print("Naming Conflict. You do not have permissions to access this organism. Either request permission from the owner, or choose a different name for your organism.")
            sys.exit(2)

        log.info("\tUpdating Organism")
        data = wa.organisms.updateOrganismInfo(
            org['id'],
            org_cn,
            args.jbrowse,
            # mandatory
            genus=args.genus,
            species=args.species,
            public=args.public
        )
        time.sleep(2)
        data = [wa.organisms.findOrganismById(org['id'])]
    else:
        # New organism
        log.info("\tAdding Organism")
        data = wa.organisms.addOrganism(
            org_cn,
            args.jbrowse,
            genus=args.genus,
            species=args.species,
            public=args.public
        )

        # Must sleep before we're ready to handle
        time.sleep(2)
        log.info("Updating permissions for %s on %s", gx_user, org_cn)
        wa.users.updateOrganismPermission(
            gx_user, org_cn,
            write=True,
            export=True,
            read=True,
        )

    data = [o for o in data if o['commonName'] == org_cn]
    print(json.dumps(data, indent=2))
