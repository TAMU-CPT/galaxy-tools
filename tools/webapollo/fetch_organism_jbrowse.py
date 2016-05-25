#!/usr/bin/env python
import os
import json
import argparse
import time
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
    org_cn = GuessOrg(args)
    org = wa.organisms.findOrganismByCn(org_cn)

    if not os.path.exists(args.target_dir):
        os.makedirs(args.target_dir)

    cmd = [
        'cp', '-R',
        org['directory'],
        os.path.join(args.target_dir, 'data')
    ]
    subprocess.check_call(cmd)
