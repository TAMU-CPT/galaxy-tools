#!/usr/bin/env python
import sys
import time
import random
import argparse
from webapollo import WAAuth, WebApolloInstance


def pwgen(length):
    chars = list("qwrtpsdfghjklzxcvbnm")
    return "".join(random.choice(chars) for _ in range(length))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sample script to add an account via web services"
    )
    WAAuth(parser)

    parser.add_argument("email", help="User Email")
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    password = pwgen(24)
    time.sleep(1)
    users = wa.users.loadUsers()

    bot_email = "bot-" + args.email.replace("@", "_") + "@cpt.tamu.edu"
    user = [u for u in users if u.username == bot_email]

    uargs = [bot_email, "BOT ACCOUNT", args.email, password]

    if len(user) == 1:
        # Update name, regen password if the user ran it again
        userObj = user[0]
        email = args.email
        q = [userObj] + uargs
        returnData = wa.users.updateUser(*q)
        sys.stdout.write("Updated User\n")
    else:
        returnData = wa.users.createUser(*uargs)
        sys.stdout.write("Created User\n")

    print("Username: %s\nPassword: %s" % (uargs[0], uargs[-1]))
    print("Return data: " + str(returnData))
