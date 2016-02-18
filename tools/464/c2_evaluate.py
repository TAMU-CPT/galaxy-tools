#!/usr/bin/env python
import argparse
import subprocess
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
from guanine import GuanineClient
g = GuanineClient()

def validate(official_data, user_data, user_email):
    od = subprocess.check_output(['md5sum', official_data])
    ud = subprocess.check_output(['md5sum', user_data])

    od = od[0:32]
    ud = ud[0:32]

    if od == ud:
        g.submit(user_email, 'C2', 1)
    else:
        g.submit(user_email, 'C2', 0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='verify against expectations')
    parser.add_argument('official_data')
    parser.add_argument('user_data')
    parser.add_argument('user_email')
    args = parser.parse_args()
    validate(**vars(args))
