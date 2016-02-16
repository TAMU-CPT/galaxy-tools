#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def guanine_pass(user_email, test_name):
    from guanine import GuanineClient
    g = GuanineClient()
    g.submit(user_email, test_name, 1.0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='verify against expectations')
    parser.add_argument('user_email', help='User email')
    parser.add_argument('test_name', help='Assessment ID')
    args = parser.parse_args()
    guanine_pass(**vars(args))
