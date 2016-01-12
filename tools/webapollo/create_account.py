#!/usr/bin/env python
import json
import random
import argparse
from webapollo import WebApolloInstance

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
    print 'Username: %s\nEmail: %s' % (args.email, password)
