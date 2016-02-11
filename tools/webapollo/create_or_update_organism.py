#!/usr/bin/env python
import json
import argparse
import time
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
    parser.add_argument('jbrowse', help='JBrowse Data Directory')
    parser.add_argument('email', help='User Email')
    parser.add_argument('--blatdb', help='BlatDB Directory')
    parser.add_argument('--genus', help='Organism Genus')
    parser.add_argument('--species', help='Organism Species')
    parser.add_argument('--public', action='store_true', help='Make organism public')

    args = parser.parse_args()


    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    gx_user = wa.users.loadUsers(email=args.email)
    if len(gx_user) == 0:
        raise Exception("Unknown user. Please register first")

    log.info("Determining if add or update required")
    try:
        org = wa.organisms.findOrganismByCn(args.cn)
    except Exception:
        org = None

    if org:
        log.info("\tUpdating Organism")
        data = wa.organisms.updateOrganismInfo(
            org['id'],
            args.cn,
            args.jbrowse,
            # mandatory
            blatdb=args.blatdb,
            genus=args.genus,
            species=args.species,
            public=args.public
        )
    else:
        # New organism
        log.info("\tAdding Organism")
        data = wa.organisms.addOrganism(
            args.cn,
            args.jbrowse,
            blatdb=args.blatdb,
            genus=args.genus,
            species=args.species,
            public=args.public
        )
        log.info("Success: %s", orgs[0]['id'])

        # Must sleep before we're ready to handle
        time.sleep(1)
        log.info("Updating permissions for %s on %s", gx_user[0], args.cn)
        wa.users.updateOrganismPermission(
            gx_user[0], args.cn,
            write=True,
            export=True,
            read=True,
        )
    data = [o for o in data if o['commonName'] == args.cn]

    print json.dumps(data, indent=2)
