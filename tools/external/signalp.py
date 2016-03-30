#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type, get_id, fetchParent
import re

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

NAME = re.compile('^# Name=(.*)\tLength=')
DATALINE = re.compile('^\s*(\d+)\s+(.)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)$')
CLEAVAGE = re.compile('^Name=(.*)\tSP=\'YES\' Cleavage site between pos. (\d+) and (\d+)')

def process(signalp):
    data = {}
    currentTarget = None
    for line in signalp:
        m = NAME.match(line)
        if m:
            currentTarget = m.groups(1)[0]
            data[currentTarget] = {
                'cleavage': False,
                'data': []
            }
            continue

        m = DATALINE.match(line)
        if m:
            position, _, c, s, y = m.groups()
            c = float(c)
            s = float(s)
            y = float(y)
            data[currentTarget]['data'].append((c, s, y))
            continue

        m = CLEAVAGE.match(line)
        if m:
            data[currentTarget]['cleavage'] = True
            data[currentTarget]['cleavageSite'] = m.groups()[1:]
            continue

    return data

def bigwig_add_header(bw_handle, identifier):
    bw_handle.write("track type=wiggle_0 name=SignalP-%s visibility=full\n" % identifier)

def bigwig_store(bw_handle, chrom, data):
    bw_handle.write("variableStep chrom=%s span=1\n" % chrom)
    for position, value in enumerate(data):
        bw_handle.write('%s %.3f\n' % (position + 1, value))


def feature_test_id(feature, **kwargs):
    return get_id(feature) in kwargs['id']

def writeGff3(data, handle, parentGff3):
    for record in GFF.parse(parentGff3):
        cdss = list(feature_lambda(
                record.features,
                feature_test_id,
                {'id': data.keys()},
                subfeatures=False
            )
        )
        record.features = []
        for cds in cdss:
            if 'note' not in cds.qualifiers:
                cds.qualifiers['note'] = []

            id = get_id(cds)
            if data[id]['cleavage']:
                cds.qualifiers['note'].append('Cleavage between %s and %s' % data[id]['cleavageSite'])
            record.features.append(fetchParent(cds))

        GFF.write([record], handle)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('parentGff3', type=file)
    parser.add_argument('signalp', type=file)
    parser.add_argument('--bw_c', default='signalp_c.wig')
    parser.add_argument('--bw_s', default='signalp_s.wig')
    parser.add_argument('--bw_y', default='signalp_y.wig')
    args = parser.parse_args()
    signalPdata = process(args.signalp)

    writeGff3(signalPdata, sys.stdout, args.parentGff3)

    bw_c = open(args.bw_c, 'w')
    bw_s = open(args.bw_s, 'w')
    bw_y = open(args.bw_y, 'w')

    bigwig_add_header(bw_c, 'c')
    bigwig_add_header(bw_s, 's')
    bigwig_add_header(bw_y, 'y')

    for sequence in signalPdata:
        # stackoverflow.com/questions/4937491/matrix-transpose-in-python
        (c, s, y) = zip(*signalPdata[sequence]['data'])

        bigwig_store(bw_c, sequence, c)
        bigwig_store(bw_s, sequence, s)
        bigwig_store(bw_y, sequence, y)

    bw_c.close()
    bw_s.close()
    bw_y.close()
