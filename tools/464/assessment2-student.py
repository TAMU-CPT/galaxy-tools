#!/usr/bin/env python
import os
import sys
import time
import uuid
import base64
import logging
import argparse
import subprocess
from webapollo import WAAuth, WebApolloInstance, AssertUser

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an account via web services')
    WAAuth(parser)
    # Need student email to share with them
    parser.add_argument('email', help='User Email')
    parser.add_argument('dataset_dir')
    args = parser.parse_args()
    # Login to Apollo
    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # We'll name our organism something easily deleteable.
    org_cn = '464-assessment-2-%s' % args.email
    # Assert that they have an account (they just need to visit now)
    gx_user = AssertUser(wa.users.loadUsers(email=args.email))
    outdir = os.path.join(args.dataset_dir, uuid.uuid4().hex)
    # Make the directory
    subprocess.check_call([
        'mkdir', '-p', outdir
    ])
    # Copy the data directory over
    subprocess.check_call([
        'tar', 'xfz', os.path.join(SCRIPT_DIR, 'test-data', 'assessment-2.tgz'),
        '--strip-components', '1', '-C', outdir
    ])
    # Now finally ready to create the organism
    data = wa.organisms.addOrganism(
        org_cn,
        outdir,
        genus="",
        species="",
        public=False
    )
    # Wait for apollo to be ready
    time.sleep(2)
    # Grant user permissions on the org
    log.info("Updating permissions for %s on %s", gx_user, org_cn)
    wa.users.updateOrganismPermission(
        gx_user, org_cn,
        write=True,
        export=True,
        read=True,
    )
    # Pulling out our results, fetch the right one from data[]
    data = [o for o in data if o['commonName'] == org_cn]

    # Generic HTML template pointing to the organism
    HTML_TPL = """
    PGh0bWw+PGhlYWQ+PHRpdGxlPkVtYmVkZGVkIEFwb2xsbyBBY2Nlc3M8L3RpdGxlPjxzdHlsZSB0
    eXBlPSJ0ZXh0L2NzcyI+Ym9keSB7e21hcmdpbjogMDt9fSBpZnJhbWUge3tib3JkZXI6IDA7d2lk
    dGg6IDEwMCU7aGVpZ2h0OiAxMDAlfX08L3N0eWxlPjwvaGVhZD48Ym9keT48aWZyYW1lIHNyYz0i
    e2Jhc2VfdXJsfS9hbm5vdGF0b3IvbG9hZExpbms/bG9jPXtjaHJvbX0mb3JnYW5pc209e29yZ0lk
    fSZ0cmFja2xpc3Q9dHJ1ZSI+PC9pZnJhbWU+PC9ib2R5PjwvaHRtbD4K
    """
    HTML_TPL = base64.b64decode(HTML_TPL.replace('\n', ''))
    # Serialize it
    sys.stdout.write(HTML_TPL.format(base_url='https://cpt.tamu.edu/apollo', chrom="", orgId=data[0]['id']))
