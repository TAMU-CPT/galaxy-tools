#!/usr/bin/env python
import json
import argparse
import time
from webapollo import WAAuth, WebApolloInstance, OrgOrGuess, GuessOrg
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)

    parser.add_argument('jbrowse', help='JBrowse Data Directory')
    parser.add_argument('email', help='User Email')
    OrgOrGuess(parser)
    parser.add_argument('--genus', help='Organism Genus')
    parser.add_argument('--species', help='Organism Species')
    parser.add_argument('--public', action='store_true', help='Make organism public')
    parser.add_argument('--secret', action='store_true', help='Make organism public')

    args = parser.parse_args()

    org_cn = GuessOrg(args)
    if isinstance(org_cn, list):
        org_cn = org_cn[0]

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    gx_user = wa.users.loadUsers(email=args.email)
    if len(gx_user) == 0:
        raise Exception("Unknown user. Please register first")

    log.info("Determining if add or update required")
    try:
        org = wa.organisms.findOrganismByCn(org_cn)
    except Exception:
        org = None

    # TODO: Check ownership
    if org:
        log.info("\tUpdating Organism")
        import subprocess
        container_id = subprocess.check_output(
            """docker ps | grep apollo | grep postgres |awk '{print $1}'""",
            shell=True
        )
        container_id = container_id.strip()

        import re
        genus = re.sub('[^A-Za-z0-9.]*', '', args.genus)
        species = re.sub('[^A-Za-z0-9.]*', '', args.species)

        SQL = """update organism set directory='%s', genus='%s', species='%s' where id='%s'""" % (
            args.jbrowse, genus, species, org['id']
        )
        CMD = """docker exec -t %s psql -U postgres -c "%s" """ % ( container_id, SQL )
        subprocess.check_output(
            CMD,
            shell=True
        )
        data = [org]
        # import sys; sys.exit()
        # data = wa.organisms.updateOrganismInfo(
            # org['id'],
            # org_cn,
            # args.jbrowse,
            # # mandatory
            # genus=args.genus,
            # species=args.species,
            # public=args.public
        # )
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
        log.info("Updating permissions for %s on %s", gx_user[0], org_cn)
        wa.users.updateOrganismPermission(
            gx_user[0], org_cn,
            write=True,
            export=True,
            read=True,
        )

    data = [o for o in data if o['commonName'] == org_cn]
    print json.dumps(data, indent=2)
