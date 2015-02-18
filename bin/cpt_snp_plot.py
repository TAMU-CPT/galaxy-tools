#!/usr/bin/env python
import logging
import random
import svgwrite
import copy
logging.basicConfig(level=logging.INFO)
import argparse
from Bio import SeqIO
import StringIO
import hashlib


def plot_snps(reference, mutated):
    reference = str(list(SeqIO.parse(reference, 'fasta'))[0].seq).upper()
    mutated = SeqIO.parse(mutated, 'fasta')

    diffs = {
    }

    for record in mutated:
        seq = str(record.seq).upper()

        if len(seq) != len(reference):
            try:
                diffs[len(seq)].append('*')
            except:
                diffs[len(seq)] = ['*']
        else:
            diff = [i for i in range(len(reference)) if reference[i] != seq[i]]
            if len(diff) > 0:
                di = diff[0]
                try:
                    diffs[di].append(seq[di])
                except:
                    diffs[di] = [seq[di]]
            else:
                pass
                # Silent mutation

    dwg = svgwrite.Drawing()
    width = 40
    global_y_offset = 20
    global_x_offset = 20
    char_width = 20
    char_height = 20

    for i in range(len(reference)/width):
        start = i * width
        end = (i+1) * width
        max_height = max([0] + [
            len(diffs[x])
            for x in diffs.keys()
            if start <= x < end])

        global_y_offset += (max_height + 4) * char_height
        # Add zero in case no muts are found

        for x_offset, char in enumerate(reference[start:end]):
            dwg.add(dwg.text(char,
                             insert=(
                                 global_x_offset + x_offset * char_width,
                                 global_y_offset
                             )))
        for x_offset, char in enumerate(reference[start:end]):
            if x_offset % 10 == 0:
                dwg.add(dwg.line(
                    start=(
                        global_x_offset + x_offset * char_width - .2 * char_width,
                        global_y_offset - 0.5 * char_height
                    ),
                    end=(
                        global_x_offset + x_offset * char_width - .2 * char_width,
                        global_y_offset + 0.5 * char_height
                    ),
                    style="stroke:black",
                ))
                dwg.add(dwg.text(x_offset + start,
                                insert=(
                                    global_x_offset + x_offset * char_width,
                                    global_y_offset + char_height
                                )))

        for j in range(start, end):
            if j in diffs:
                for k, char in enumerate(sorted(diffs[j])):
                    dwg.add(dwg.text(char,
                                     insert=(
                                         global_x_offset + (j - start) * char_width,
                                         global_y_offset - (1 + k) * char_height
                                     )))


    print dwg.tostring()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='simple plot of generated snps against a reference sequence')
    parser.add_argument('reference', type=file, help='Reference sequence (only 1 fasta sequence allowed)')
    parser.add_argument('mutated', type=file, help='Fasta file of SNPs')

    args = parser.parse_args()
    plot_snps(**vars(args))
