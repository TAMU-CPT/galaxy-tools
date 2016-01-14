#!/usr/bin/env python
import random
import argparse
import time
from webapollo import WebApolloInstance, GroupObj

def pwgen(length):
    chars = list('qwrtpsdfghjklzxcvbnm')
    return ''.join(random.choice(chars) for _ in range(length))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an account via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Admin Username')
    parser.add_argument('password', help='WA Admin Password')

    parser.add_argument('email', help='User Email')
    parser.add_argument('--first', help='First Name', default='J')
    parser.add_argument('--last', help='Last Name', default='Aggie')
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    password = pwgen(12)
    wa.users.createUser(args.email, args.first, args.last, password, role='user')
    time.sleep(1)
    user = [u for u in wa.users.loadUsers()
            if u.username == args.email][0]

    bich464 = GroupObj(name="bich464-2016-spring")

    # Update name, regen password if the user ran it again
    wa.users.updateUser(user, args.email, args.first, args.last, password)
    # Add to bich464 group
    wa.users.addUserToGroup(bich464, user)

    print 'Username: %s\nPassword: %s' % (args.email, password)
