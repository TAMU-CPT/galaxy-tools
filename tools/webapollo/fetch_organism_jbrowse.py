#!/usr/bin/env python
import os
import sys
import time
import argparse
from webapollo import WAAuth, WebApolloInstance, GuessOrg, OrgOrGuess
import logging
import subprocess
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)
    OrgOrGuess(parser)
    parser.add_argument('target_dir', help='Target directory')

    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    org_cn = GuessOrg(args, wa)
    if isinstance(org_cn, list):
        org_cn = org_cn[0]
    org = wa.organisms.findOrganismByCn(org_cn)

    if not os.path.exists(args.target_dir):
        os.makedirs(args.target_dir)

    if not os.path.exists(os.path.join(org['directory'], 'seq')):
        sys.stderr.write("Missing seq directory BEFORE copy")
        sys.exit(1)

    cmd = [
        'rsync', '-avr',
        org['directory'].rstrip('/') + '/',
        os.path.join(args.target_dir, 'data', '')
    ]
    # We run this OBSESSIVELY because my org had a hiccup where the origin
    # (silent) cp -R failed at one point. This caused MANY HEADACHES.
    #
    # Our response is to run this 3 times (in case the issue is temporary),
    # with delays in between. And ensure that we have the correct number of
    # files / folders before and after.
    sys.stderr.write(' '.join(cmd))
    sys.stderr.write('\n')
    sys.stderr.write(subprocess.check_output(cmd))
    time.sleep(5)
    sys.stderr.write('\n')
    sys.stderr.write('\n')
    sys.stderr.write(' '.join(cmd))
    sys.stderr.write('\n')
    sys.stderr.write(subprocess.check_output(cmd))
    time.sleep(5)
    sys.stderr.write('\n')
    sys.stderr.write('\n')
    sys.stderr.write(' '.join(cmd))
    sys.stderr.write('\n')
    sys.stderr.write(subprocess.check_output(cmd))

    if not os.path.exists(os.path.join(args.target_dir, 'data', 'seq')):
        sys.stderr.write("ERROR: Missing seq directory AFTER copy")
        sys.exit(1)
