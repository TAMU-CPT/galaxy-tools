#!/usr/bin/env python
import os
import json
import argparse
import time
from webapollo import WebApolloInstance
import logging
import subprocess
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Admin Username')
    parser.add_argument('password', help='WA Admin Password')
    parser.add_argument('cn', help='Organism Common Name')
    parser.add_argument('email', help='User Email')
    parser.add_argument('target_dir', help='Target directory')

    args = parser.parse_args()


    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    gx_user = wa.users.loadUsers(email=args.email)
    if len(gx_user) == 0:
        raise Exception("Unknown user. Please register first")
    org = wa.organisms.findOrganismByCn(args.cn)


    if not os.path.exists(args.target_dir):
        os.makedirs(args.target_dir)

    subprocess.check_call([
        'cp', '-R',
        os.path.join(org['directory'], '.'),
        args.target_dir
    ])
