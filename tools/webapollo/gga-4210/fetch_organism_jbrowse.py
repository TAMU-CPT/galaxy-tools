#!/usr/bin/env python
from __future__ import print_function

import argparse
import filecmp
import logging
import os
import subprocess
import sys
import time

from apollo import accessible_organisms
from apollo.util import GuessOrg, OrgOrGuess

from arrow.apollo import get_apollo_instance

from webapollo import UserObj, handle_credentials

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def are_dir_trees_equal(dir1, dir2):
    """
    Compare two directories recursively. Files in each directory are
    assumed to be equal if their names and contents are equal.

    @param dir1: First directory path
    @param dir2: Second directory path

    @return: True if the directory trees are the same and
        there were no errors while accessing the directories or files,
        False otherwise.

    # http://stackoverflow.com/questions/4187564/recursive-dircmp-compare-two-directories-to-ensure-they-have-the-same-files-and/6681395#6681395
    """

    dirs_cmp = filecmp.dircmp(dir1, dir2)
    if len(dirs_cmp.left_only) > 0 or len(dirs_cmp.right_only) > 0 or \
            len(dirs_cmp.funny_files) > 0:
        print(('LEFT', dirs_cmp.left_only))
        print(('RIGHT', dirs_cmp.right_only))
        print(('FUNNY', dirs_cmp.funny_files))
        return False
    (_, mismatch, errors) = filecmp.cmpfiles(
        dir1, dir2, dirs_cmp.common_files, shallow=False)
    if len(mismatch) > 0 or len(errors) > 0:
        print(mismatch)
        print(errors)
        return False
    for common_dir in dirs_cmp.common_dirs:
        new_dir1 = os.path.join(dir1, common_dir)
        new_dir2 = os.path.join(dir2, common_dir)
        if not are_dir_trees_equal(new_dir1, new_dir2):
            return False
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to fetch organism data dir from apollo')
    OrgOrGuess(parser)
    parser.add_argument('target_dir', help='Target directory')
    parser.add_argument('email', help='User Email')

    args = parser.parse_args()

    wa = get_apollo_instance()
    # User must have an account
    org_cn = GuessOrg(args, wa)
    if isinstance(org_cn, list):
        org_cn = org_cn[0]

    # User must have an account, if not, create it
    gx_user = UserObj(**wa.users._assert_or_create_user(args.email))
    handle_credentials(gx_user)

    all_orgs = wa.organisms.get_organisms()
    if 'error' in all_orgs:
        all_orgs = []
    all_orgs = [org['commonName'] for org in all_orgs]
    if org_cn not in all_orgs:
        raise Exception("Could not find organism %s" % org_cn)

    # User must have READ access
    orgs = accessible_organisms(gx_user, [org_cn], 'READ')
    if not orgs:
        raise Exception("You do not have write permission on this organism")
    org = wa.organisms.show_organism(org_cn)

    if not os.path.exists(args.target_dir):
        os.makedirs(args.target_dir)

    if not os.path.exists(os.path.join(org['directory'], 'seq')):
        sys.stderr.write("Missing seq directory BEFORE copy")
        sys.exit(1)

    cmd = [
        'rsync', '-avr',
        org['directory'].rstrip('/') + '/',
        os.path.join(args.target_dir, 'data', '')
    ]
    # We run this OBSESSIVELY because my org had a hiccup where the origin
    # (silent) cp -R failed at one point. This caused MANY HEADACHES.
    #
    # Our response is to run this 3 times (in case the issue is temporary),
    # with delays in between. And ensure that we have the correct number of
    # files / folders before and after.
    sys.stderr.write(' '.join(cmd))
    sys.stderr.write('\n')
    sys.stderr.write(subprocess.check_output(cmd).decode(sys.stderr.encoding))
    if not are_dir_trees_equal(
        os.path.join(org['directory'].rstrip('/')),
        os.path.join(args.target_dir, 'data')
    ):
        # Not good
        time.sleep(5)
        sys.stderr.write('\n')
        sys.stderr.write(' '.join(cmd))
        sys.stderr.write('\n')
        sys.stderr.write(subprocess.check_output(cmd).decode(sys.stderr.encoding))
        if not are_dir_trees_equal(
            os.path.join(org['directory'].rstrip('/'), 'data'),
            os.path.join(args.target_dir, 'data')
        ):
            time.sleep(5)
            sys.stderr.write('\n')
            sys.stderr.write(' '.join(cmd))
            sys.stderr.write('\n')
            sys.stderr.write(subprocess.check_output(cmd).decode(sys.stderr.encoding))
            if not are_dir_trees_equal(
                os.path.join(org['directory'].rstrip('/'), 'data'),
                os.path.join(args.target_dir, 'data')
            ):
                sys.stderr.write('FAILED THREE TIMES TO COPY. SOMETHING IS WRONG WRONG WRONG.')
                sys.exit(2)
