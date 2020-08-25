#!/usr/bin/env python
from __future__ import print_function

import argparse
import json

from apollo import accessible_organisms

from arrow.apollo import get_apollo_instance

from webapollo import UserObj, handle_credentials

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='List all organisms available in an Apollo instance')
    parser.add_argument('email', help='User Email')
    args = parser.parse_args()

    wa = get_apollo_instance()

    gx_user = UserObj(**wa.users._assert_or_create_user(args.email))
    handle_credentials(gx_user)

    all_orgs = wa.organisms.get_organisms()

    orgs = accessible_organisms(gx_user, all_orgs)

    print(json.dumps(orgs, indent=2))
