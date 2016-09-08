#!/usr/bin/env python
from BCBio import GFF
import argparse
import sys

def parse_gff(id_start_end, gff3):
    ids = {}
    for line in id_start_end:
        l = line.split()
        if l[0] in ids:
            ids[l[0]].append((int(l[1]), int(l[2])))
        else:
            ids[l[0]] = [(int(l[1]), int(l[2]))]

    for rec in GFF.parse(gff3):
        locs = ids[rec.id]
        feats = []
        for feat in rec.features:
            f_loc = (feat.location.start, feat.location.end)
            for loc in locs:
                if (min(f_loc) <= max(loc)) and (max(f_loc) >= min(loc)):
                    feats.append(feat)

        rec.features = feats
        GFF.write([rec], sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='features from a GFF3 file based on location')
    parser.add_argument('id_start_end', type=file, help='table with record ids, start locations, and end locations')
    parser.add_argument('gff3', type=file, help='GFF3 annotations')
    args = parser.parse_args()
    parse_gff(**vars(args))
