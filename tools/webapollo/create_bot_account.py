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
    parser.add_argument('--first', help='First Name', default='Jane')
    parser.add_argument('--last', help='Last Name', default='Aggie')
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    password = pwgen(24)
    time.sleep(1)
    users = wa.users.loadUsers()
    user = [u for u in users
            if u.username == args.email]

    if len(user) == 1:
        # Update name, regen password if the user ran it again
        userObj = user[0]
        returnData = wa.users.updateUser(userObj, args.email, args.first, args.last, password)
        print 'Updated User\nUsername: %s\nPassword: %s' % (args.email, password)
    else:
        returnData = wa.users.createUser(args.email, args.first, args.last, password, role='user')
        print 'Created User\nUsername: %s\nPassword: %s' % (args.email, password)

    print "Return data: " + str(returnData)
