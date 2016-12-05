#!/usr/bin/env python
import numpy
import argparse
from phantasm import Utils, Clustering

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Columnar Data Rescale and Partitioning')
    parser.add_argument('tabular_data', type=argparse.FileType("r"), help='Tabular Dataset')

    choices = ['MeanShift', 'user', 'simple']
    parser.add_argument('partition_type', choices=choices, help='Partitioning Method')

    parser.add_argument('--user_breakpoints', type=float, nargs='*', help='User specified breakpoints')
    parser.add_argument('--percentage_breakpoints', type=int, nargs='1',
                        help='Percentage based breakpoints into N bins')

    parser.add_argument('--version', action='version', version='0.2')
    args = parser.parse_args()

    data = Utils.load_data(args.tabular_data)
    (id_col, data_col) = (data['id'], data['data'])

    if args.partition_type == 'MeanShift':
        data_col = Clustering.meanshift(data_col)
    elif args.partition_type == 'user':
        data_col = Clustering.user_break(data_col, breaks=args.user_breakpoints)
    elif args.partition_type == 'simple':
        dataset_min = numpy.min(data_col)
        dataset_max = numpy.max(data_col)
        break_size = (dataset_max - dataset_min) / args.percentage_breakpoints
        break_points = range(dataset_min, dataset_max + break_size, break_size)
        print dataset_min
        print dataset_max
        print break_points
        data_col = Clustering.user_break(data_col, breaks=break_points, exact=True)
    else:
        raise NotImplementedError()

    print "# ID\tData"
    for (id_val, data_val) in zip(id_col, data_col):
        print "%s\t%s" % (id_val, data_val)
