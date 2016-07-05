#!/usr/bin/env python
import argparse
import Levenshtein
from phantasm import Cassettes, Utils
import logging
logging.basicConfig(level=logging.INFO)


def compare_values(a, b, method, undef_value):
    if method not in ('phantasm_cids', 'texteq', 'levenshtein'):
        # Methods listed above use strings, so we can ignore those.
        #
        # If we can't convert to float, then we just exit early and return our
        # undefined value
        try:
            a = float(a)
        except:
            return undef_value
        try:
            b = float(b)
        except:
            return undef_value

    if method == 'add':
        return a + b
    elif method == 'sub':
        return a - b
    elif method == 'mult':
        return a * b
    elif method == 'div':
        if b == 0:
            return undef_value
        else:
            return a / b
    elif method == 'dist':
        return abs(a - b)
    elif method == 'pdiff':
        return abs(a - b) / (a + b)
    elif method == 'bit_diff':
        a = abs(int(a))
        b = abs(int(b))
        return bin(a ^ b).count("1")
    elif method == 'phantasm_cids':
        chunked = Cassettes.revcomrot(a)
        # Compare each for levenshtein distance, return the minimum
        scores = [abs(Levenshtein.distance(x, b)) for x in chunked]
        return min(scores)
    elif method == 'numeq':
        return 1 if a == b else 0
    elif method == 'texteq':
        a = a.strip()
        b = b.strip()

        if a is None or b is None:
            return 0
        elif len(a) == 0 or len(b) == 0:
            return 0
        else:
            return 1 if a == b else 0
    elif method == 'levenshtein':
        return Levenshtein.distance(a, b)
    else:
        raise NotImplementedError("Comparison method %s not available" % method)
    return None


def compare_files(file_a, file_b, comparison_method, undef_value, **kwargs):
    (header_a, data_a) = Utils.load_data_with_headers(file_a)
    (header_b, data_b) = Utils.load_data_with_headers(file_b)

    header = ['# ID_A', 'ID_B', 'Value']
    yield header

    for row_a in data_a:
        for row_b in data_b:
            yield (row_a[0], row_b[0],
                   compare_values(row_a[1], row_b[1], comparison_method, undef_value))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate comparisons between two PHAnTASM datasets')
    parser.add_argument('file_a', type=file, help='First Tabular Dataset')
    parser.add_argument('file_b', type=file, help='Second Tabular Dataset')
    choices = ['pdiff', 'bit_diff', 'add', 'sub', 'mult', 'div', 'dist', 'phantasm_cids',
               'levenshtein', 'texteq', 'numeq']
    parser.add_argument('comparison_method', choices=choices,
                        help='Method for comparison')
    parser.add_argument('undef_value', nargs='?', type=float, default=0,
                        help='Undefined value. For operations involving division, ' +
                        'what should undefined results be set to? (e.g. 3/0 = ?).')
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()

    for row in compare_files(**vars(args)):
        print '\t'.join(map(str, row))
