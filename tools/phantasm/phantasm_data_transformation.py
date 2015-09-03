#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO)
from phantasm import Utils, Transforms

if __name__ == '__main__':
    transformations = {
        'none': 'x',
        'neg': '-x',
        'log': 'log(x)',
        'ln': 'ln(x)',
        'inv': '1/x',
        'abs': 'abs',
        'exp': 'exp',
    }

    parser = argparse.ArgumentParser(description='Transform data')
    parser.add_argument('tabular_data', type=file, help='Tabular Dataset')
    parser.add_argument('--operation', choices=transformations.keys(),
                        nargs='+', help='Transformation to apply to the data')
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()

    data = Utils.load_data(args.tabular_data)
    (id_col, data_col) = (data['id'], data['data'])

    for operation in args.operation:
        data_col = Transforms.apply_transform(data_col, operation)

    print "# ID\tData"
    for (id_val, data_val) in zip(id_col, data_col):
        print "%s\t%s" % (id_val, data_val)
