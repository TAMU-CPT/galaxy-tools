#!/usr/bin/env python
import os
import json
import argparse
import time
from webapollo import WAAuth, WebApolloInstance
import logging
import subprocess
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)

    parser.add_argument('cn', help='Organism Common Name')
    parser.add_argument('target_dir', help='Target directory')

    args = parser.parse_args()


    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    org = wa.organisms.findOrganismByCn(args.cn)

    if not os.path.exists(args.target_dir):
        os.makedirs(args.target_dir)

    cmd = [
        'cp', '-R',
        org['directory'],
        os.path.join(args.target_dir, 'data')
    ]
    subprocess.check_call(cmd)
