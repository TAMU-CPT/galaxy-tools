#!/usr/bin/env python
import os
import shutil
import argparse
import tempfile
import subprocess
# File magic numbers to determine archive type
import magic
# Util to extract archives
from archive import Archive


def main():
    #Parse Command Line
    parser = argparse.ArgumentParser(
        description='Compile an archive of circos data and config')
    parser.add_argument('--archive', dest='archive',
                        help='Use this archive for the image generation',
                        required=True)
    parser.add_argument('--png', dest='png', help='png output',
                        default='circos.png')
    parser.add_argument('--svg', dest='svg', help='svg output',
                        default='circos.svg')

    options = parser.parse_args()
    tmpdir = extract_archive(options.archive)
    circos_conf = find_conf_location(tmpdir)
    run_circos(options, circos_conf, tmpdir)


def extract_archive(archive_location):
    """Extract an archive (magically guessing filetype) into a temporary
    directory and return the temp dir's location"""
    tmpdir = tempfile.mkdtemp()
    arc = Archive(archive_location)
    arc.extract(tmpdir)
    print "Extracted to %s" % (tmpdir,)
    return tmpdir


def find_conf_location(tmpdir):
    """Given a temporary directory where an archive of circos conf and data has
    been extracted, attempt to locate the main circos.conf file"""
    import fnmatch
    matches = []
    for root, dirnames, filenames in os.walk(tmpdir):
        for filename in fnmatch.filter(filenames, "circos.conf"):
            matches.append(os.path.join(root, filename))
    if len(matches) == 1:
        return matches[0]
    elif len(matches) > 1:
        raise Exception("Found multiple circos.conf files")
    else:
        raise Exception("Could not find main circos.conf file")


def run_circos(options, conf_loc, tmpdir):
    """Run circos on a given conf_location

    tmpdir allows for cleanup afterwards, and options must have options.svg and
    options.png defined"""
    args = ['/opt/circos-0.66/bin/circos', '-conf', conf_loc]

    proc = subprocess.Popen(" ".join(args), shell=True)
    proc.communicate()

    # generate file locations
    png_from = os.path.join(os.path.dirname(conf_loc), 'circos.png')
    svg_from = os.path.join(os.path.dirname(conf_loc), 'circos.svg')

    # Copy files over
    if os.path.exists(png_from):
        shutil.move(png_from, options.png)
    if os.path.exists(svg_from):
        shutil.move(svg_from, options.svg)
    # Cleanup
    shutil.rmtree(tmpdir)


if __name__ == '__main__':
    main()
