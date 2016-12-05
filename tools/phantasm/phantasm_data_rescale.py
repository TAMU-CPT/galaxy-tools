#!/usr/bin/env python
import argparse
import logging
from phantasm import Utils, Transforms
logging.basicConfig(level=logging.INFO)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rescale data')
    parser.add_argument('tabular_data', type=argparse.FileType("r"), help='Tabular Dataset')
    parser.add_argument('--min', type=float, help='New Minimum', default=0.0)
    parser.add_argument('--max', type=float, help='New Maximum', default=1.0)

    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()

    data = Utils.load_data(args.tabular_data)
    (id_col, data_col) = (data['id'], data['data'])
    (squash, offset) = Transforms.calculate_reshape(data_col, args.min, args.max)
    data_col = Transforms.apply_reshape(data_col, squash, offset)

    print "# ID\tData"
    for (id_val, data_val) in zip(id_col, data_col):
        print "%s\t%s" % (id_val, data_val)
