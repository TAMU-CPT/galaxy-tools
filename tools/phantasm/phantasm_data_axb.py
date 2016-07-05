#!/usr/bin/env python
import argparse
import logging
import numpy
logging.basicConfig(level=logging.INFO)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Transform data')
    parser.add_argument('tabular_data', type=file, help='Tabular Dataset')
    parser.add_argument('m', type=float, help='m of mx+b')
    parser.add_argument('b', type=float, help='b of mx+b')
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()

    data = numpy.genfromtxt(args.tabular_data, delimiter='\t', dtype=None, names=('id1', 'id2', 'data'))
    (id_col, data_col) = (data['id'], data['data'])

    data_col = (data_col * args.m) + args.b

    print "# ID\tData"
    for (id_val, data_val) in zip(id_col, data_col):
        print "%s\t%s" % (id_val, data_val)
