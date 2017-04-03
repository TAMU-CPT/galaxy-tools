#!/usr/bin/env python
import os
import sys
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
        org['directory'].rstrip('/'),
        os.path.join(args.target_dir, 'data')
    ]
    sys.stderr.write(' '.join(cmd))
    sys.stderr.write(subprocess.check_output(cmd))
    sys.stderr.write(' '.join(cmd))
    sys.stderr.write(subprocess.check_output(cmd))
    sys.stderr.write(' '.join(cmd))
    sys.stderr.write(subprocess.check_output(cmd))

    if not os.path.exists(os.path.join(args.target_dir, 'data', 'seq')):
        sys.stderr.write("Missing seq directory AFTER copy")
        sys.exit(1)
