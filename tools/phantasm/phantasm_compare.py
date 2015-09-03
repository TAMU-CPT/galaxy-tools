#!/usr/bin/env python
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PHAnTASM Comparison Mapper')
    parser.add_argument('--files', type=file, nargs="+", help='Input Two-Way Comparison')
    parser.add_argument('--weights', type=float, nargs="+", help="Metric Weighting")
    parser.add_argument('--version', action='version', version='0.2')
    args = parser.parse_args()

    if len(args.files) != len(args.weights):
        raise Exception("Must specifiy same number of files and metric_weights")

    result_map = {}
    weight_sums = sum(args.weights)

    for f, m in zip(args.files, args.weights):
        data = [x.strip().split('\t') for x in f.readlines() if not
                x.startswith('#')]
        for (a, b, score) in data:
            if a not in result_map:
                result_map[a] = {}
            if b not in result_map[a]:
                result_map[a][b] = 0

            result_map[a][b] += (float(score) * m/weight_sums)

    keys = sorted(result_map.keys())
    # Header
    print '\t'.join([''] + keys)
    for f in keys:
        current_row = [f]
        for t in keys:
            current_row.extend([result_map[f][t]])
        print '\t'.join(current_row)
