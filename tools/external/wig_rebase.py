#!/usr/bin/env python
import re
import sys
import logging
import argparse
import numpy
from gff3 import feature_lambda, feature_test_true
from BCBio import GFF
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def __update_feature_location(pos, parent, protein2dna):
    if protein2dna:
        pos *= 3
        # Move back so the "span" is correct.

    # print(start, end, parent.location.start, parent.location.end)
    if parent.strand >= 0:
        new_pos = parent.start + pos
    else:
        new_pos = parent.end - pos

    # print(start, end, ns, ne, st)

    # Don't let start/stops be less than zero. It's technically valid for them
    # to be (at least in the model I'm working with) but it causes numerous
    # issues.
    #
    # Instead, we'll replace with %3 to try and keep it in the same reading
    # frame that it should be in.
    if new_pos < 0:
        new_pos %= 3
    return new_pos


def getGff3Locations(parent, map_by='ID'):
    featureLocations = {}
    recs = GFF.parse(parent)
    # Only parse first.
    rec = next(recs)
    # Get all the feature locations in this genome
    for feature in feature_lambda(rec.features, feature_test_true, {}):
        id = feature.qualifiers.get(map_by, [feature.id])[0]
        featureLocations[id] = feature.location
    return rec, featureLocations


def rebase_wig(parent, wigData, protein2dna=False, map_by='ID'):
    rec, locations = getGff3Locations(parent, map_by=map_by)
    current_id = None
    current_ft = None
    # We have to store in a giant array so we can overwrite safely and don't
    # emit multiple values. UCSC tools are sucky.
    values = numpy.zeros(500000)
    for line in wigData:
        if line.startswith('track'):
            # pass through
            sys.stdout.write(line)
            sys.stdout.write('variableStep chrom=%s span=%s\n' % (rec.id, 1 if not protein2dna else 3))
            continue
        if line.startswith('variableStep'):
            # No passthrough
            current_id = re.findall('chrom=([^ ]+)', line)[0]
            current_ft = locations[current_id]
        else:
            (pos, val) = line.strip().split()
            pos = int(pos)
            val = float(val)

            npos = __update_feature_location(pos, current_ft, protein2dna=protein2dna)
            values[npos] = val

    for i in range(len(values)):
        if values[i] != 0:
            sys.stdout.write('%s %s\n' % (i, values[i]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase wig data against parent locations')
    parser.add_argument('parent', type=argparse.FileType("r"))
    parser.add_argument('wigData', type=argparse.FileType("r"))
    parser.add_argument('--protein2dna', action='store_true',
                        help='Map protein translated results to original DNA data')
    parser.add_argument('--map_by', help='Map by key', default='ID')
    args = parser.parse_args()
    rebase_wig(**vars(args))
