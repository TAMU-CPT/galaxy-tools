#!/usr/bin/env python
import random
import argparse
import time
from webapollo import WAAuth, WebApolloInstance


def pwgen(length):
    chars = list('qwrtpsdfghjklzxcvbnm')
    return ''.join(random.choice(chars) for _ in range(length))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an account via web services')
    WAAuth(parser)

    parser.add_argument('email', help='User Email')
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    password = pwgen(24)
    # ??
    # time.sleep(1)
    users = wa.users.loadUsers()
    user = [u for u in users
            if u.username == args.email]

    # Create user if does not exist, updating existing users is currently unhappy
    if len(user) == 0:
        returnData = wa.users.createUser(args.email, 'REMOTE', 'USER', password, role='user')

    time.sleep(1)
    users = wa.users.loadUsers()
    user = [u for u in users
            if u.username == args.email]

    print(wa.users.updateOrganismPermission(
        user[0],
        '464_2017_assessment1',
        write=True,
        read=True,
        export=True
    ))
