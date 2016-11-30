#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def check_pw(password, user_email):
    if password == 'Y6h2qhVTwWmYjFgv2V':
        pass
    else:
        import sys
        sys.exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='verify against expectations')
    parser.add_argument('user_email', help='User email')
    parser.add_argument('password', help='Password')
    args = parser.parse_args()
    check_pw(**vars(args))
