#!/usr/bin/env python
import os
import time
import base64
import logging
import argparse
from webapollo import WAAuth, WebApolloInstance, AssertUser

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an account via web services')
    WAAuth(parser)
    parser.add_argument('email', help='User Email')

    args = parser.parse_args()
    wa = WebApolloInstance(args.apollo, args.username, args.password)

    org_cn = '464-assessment-2-%s' % args.email
    gx_user = AssertUser(wa.users.loadUsers(email=args.email))

    data = wa.organisms.addOrganism(
        org_cn,
        args.jbrowse,
        genus=args.genus,
        species=args.species,
        public=args.public
    )

    time.sleep(2)
    log.info("Updating permissions for %s on %s", gx_user, org_cn)
    wa.users.updateOrganismPermission(
        gx_user, org_cn,
        write=True,
        export=True,
        read=True,
    )

    data = [o for o in data if o['commonName'] == org_cn]


    HTML_TPL = """
    PGh0bWw+PGhlYWQ+PHRpdGxlPkVtYmVkZGVkIEFwb2xsbyBBY2Nlc3M8L3RpdGxlPjxzdHlsZSB0
    eXBlPSJ0ZXh0L2NzcyI+Ym9keSB7e21hcmdpbjogMDt9fSBpZnJhbWUge3tib3JkZXI6IDA7d2lk
    dGg6IDEwMCU7aGVpZ2h0OiAxMDAlfX08L3N0eWxlPjwvaGVhZD48Ym9keT48aWZyYW1lIHNyYz0i
    e2Jhc2VfdXJsfS9hbm5vdGF0b3IvbG9hZExpbms/bG9jPXtjaHJvbX0mb3JnYW5pc209e29yZ0lk
    fSZ0cmFja2xpc3Q9dHJ1ZSI+PC9pZnJhbWU+PC9ib2R5PjwvaHRtbD4K
    """
    HTML_TPL = base64.b64decode(HTML_TPL.replace('\n', ''))

    print HTML_TPL.format(base_url='https://cpt.tamu.edu/apollo', chrom="", orgId=data[0]['id'])
