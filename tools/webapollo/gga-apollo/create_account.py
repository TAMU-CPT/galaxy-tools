#!/usr/bin/env python
from __future__ import print_function

import argparse
import time

from arrow.apollo import get_apollo_instance


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an account via web services')
    parser.add_argument('email', help='User Email')
    parser.add_argument('--first', help='First Name', default='Jane')
    parser.add_argument('--last', help='Last Name', default='Aggie')
    args = parser.parse_args()

    wa = get_apollo_instance()

    password = wa.users._password_generator(12)
    time.sleep(1)
    users = wa.users.get_users()
    user = [u for u in users
            if u['username'] == args.email]

    if len(user) == 1:
        # Update name, regen password if the user ran it again
        returnData = wa.users.update_user(args.email, args.first, args.last, password)
        print('Updated User\nUsername: %s\nPassword: %s' % (args.email, password))
    else:
        returnData = wa.users.create_user(args.email, args.first, args.last, password, role='user')
        print('Created User\nUsername: %s\nPassword: %s' % (args.email, password))

    print("Return data: " + str(returnData))
