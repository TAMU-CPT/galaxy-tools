#!/usr/bin/env python
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
import itertools
import pygal
from pygal.style import LightSolarizedStyle


def tntable(table=11):
    codons = itertools.product('ACTG', repeat=3)
    _table = {}
    translation = {}
    for x in [''.join(y) for y in codons]:
        aa = str(Seq(x, IUPAC.IUPACUnambiguousDNA()).translate(table=table))
        translation[x] = aa

        try:
            _table[aa][x] = {}
        except:
            _table[aa] = {x: {}}
    return _table, translation


def main(reference, comparison):
    ref_data = np.genfromtxt(reference, dtype=None, delimiter='\t', names=True)
    comp_data = np.genfromtxt(comparison, dtype=None, delimiter='\t', names=True)
    table, translation = tntable()

    for codon in ref_data:
        # AAA \t N
        # N -> AAA = 32
        nt = codon[0]
        aa = translation[nt]
        table[aa][nt]['ref'] = codon[1]

    for codon in comp_data:
        nt = codon[0]
        aa = translation[nt]
        table[aa][nt]['comp'] = codon[1]

    keys = []
    ref_flat = []
    comp_flat = []
    for x in sorted(table.keys()):
        for y in sorted(table[x].keys()):
            keys.append('%s (%s)' % (x, y))
            ref_val = table.get(x, {}).get(y, {}).get('ref', None)
            comp_val = table.get(x, {}).get(y, {}).get('comp', None)
            ref_flat.append({'value': ref_val, 'label': str(ref_val)})
            comp_flat.append({'value': comp_val, 'label': str(comp_val)})

    bar_chart = pygal.Bar(width=1500, height=400, x_label_rotation=90, style=LightSolarizedStyle, delete=False)
    bar_chart.title = 'Codon usage comparison'
    bar_chart.x_labels = keys
    bar_chart.add('Reference', ref_flat)
    bar_chart.add('Comparison', comp_flat)
    print bar_chart.render()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot comparison of codon usage', epilog="")
    parser.add_argument('reference', type=file, help='Reference Data')
    parser.add_argument('comparison', type=file, help='Data to compare')
    args = parser.parse_args()
    main(**vars(args))
