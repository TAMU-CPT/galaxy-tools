#!/usr/bin/env python
import argparse
from webapollo import WebApolloInstance, AssertAdmin, AssertUser


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sample script to add an attribute to a feature via web services"
    )
    parser.add_argument("apollo", help="Complete Apollo URL")
    parser.add_argument("username", help="WA Admin Username")
    parser.add_argument("password", help="WA Admin Password")

    parser.add_argument("email", help="Your email")
    args = parser.parse_args()
    wa = WebApolloInstance(args.apollo, args.username, args.password)
    gx_user = AssertAdmin(AssertUser(wa.users.loadUsers(email=args.email)))

    for user in wa.users.loadUsers():
        print(
            "%s\t%s\t%s\t%s"
            % (user.userId, user.firstName, user.lastName, user.username)
        )
