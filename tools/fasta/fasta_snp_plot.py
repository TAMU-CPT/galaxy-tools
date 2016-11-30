#!/usr/bin/env python
import logging
import svgwrite
import argparse
from Bio import SeqIO
logging.basicConfig(level=logging.INFO)


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
    width = 120
    global_y_offset = 20
    global_x_offset = 20
    char_width = 10
    char_height = 15

    color_ref = {
        'r': ('A91a00', 'CA290C'),
        'g': ('045069', '0C627E'),
        'o': ('A97700', 'ca920c'),
    }

    colors = {
        # Charged:
        'R': color_ref['r'][0], 'K': color_ref['r'][0], 'D': color_ref['r'][0], 'E': color_ref['r'][0],
        # Polar (may participate in hydrogen bonds):
        'Q': color_ref['g'][0], 'N': color_ref['g'][0], 'H': color_ref['g'][0], 'S': color_ref['g'][0],
        'T': color_ref['g'][0], 'Y': color_ref['g'][0], 'C': color_ref['g'][0], 'M': color_ref['g'][0],
        'W': color_ref['g'][0],
        # Hydrophobic (normally buried inside the protein core):
        'A': color_ref['o'][0], 'I': color_ref['o'][0], 'L': color_ref['o'][0], 'F': color_ref['o'][0],
        'V': color_ref['o'][0], 'P': color_ref['o'][0], 'G': color_ref['o'][0],
    }
    colors_muted = {
        # Charged:
        'R': color_ref['r'][1], 'K': color_ref['r'][1], 'D': color_ref['r'][1], 'E': color_ref['r'][1],
        # Polar (may participate in hydrogen bonds):
        'Q': color_ref['g'][1], 'N': color_ref['g'][1], 'H': color_ref['g'][1], 'S': color_ref['g'][1],
        'T': color_ref['g'][1], 'Y': color_ref['g'][1], 'C': color_ref['g'][1], 'M': color_ref['g'][1],
        'W': color_ref['g'][1],
        # Hydrophobic (normally buried inside the protein core):
        'A': color_ref['o'][1], 'I': color_ref['o'][1], 'L': color_ref['o'][1], 'F': color_ref['o'][1],
        'V': color_ref['o'][1], 'P': color_ref['o'][1], 'G': color_ref['o'][1],
    }
    yrange = range(len(reference) / width)
    if len(yrange) == 0:
        yrange = [0]

    for i in yrange:
        start = i * width
        end = (i + 1) * width
        max_height = max([0] + [
            len(diffs[x])
            for x in diffs.keys()
            if start <= x < end])

        global_y_offset += (max_height + 4) * char_height
        # Add zero in case no muts are found

        # Reference sequence
        for x_offset, char in enumerate(reference[start:end]):
            dwg.add(dwg.text(char,
                             insert=(
                                 global_x_offset + x_offset * char_width,
                                 global_y_offset
                             ),
                             style="font-family:monospace;fill:#%s" % colors.get(char, 'black'),
                             ))

        # Scale
        for x_offset, char in enumerate(reference[start:end]):
            if x_offset % 10 == 0:
                dwg.add(dwg.line(
                    start=(
                        global_x_offset + x_offset * char_width - .1 * char_width,
                        global_y_offset
                    ),
                    end=(
                        global_x_offset + x_offset * char_width - .1 * char_width,
                        global_y_offset + char_height
                    ),
                    style="stroke:black",
                ))
                dwg.add(dwg.text(
                    x_offset + start,
                    insert=(
                        global_x_offset + x_offset * char_width,
                        global_y_offset + char_height
                    ),
                    style="font-family:monospace"))

        # SNPs
        for j in range(start, end):
            if j in diffs:
                for k, char in enumerate(sorted(diffs[j])):
                    dwg.add(dwg.text(char,
                                     insert=(
                                         global_x_offset + (j - start) * char_width,
                                         global_y_offset - (1 + k) * char_height
                                     ),
                                     style="font-family:monospace; fill:#%s" % colors_muted.get(char, 'black'),
                                     ))

    print dwg.tostring()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='simple plot of generated snps against a reference sequence')
    parser.add_argument('reference', type=file, help='Reference sequence (only 1 fasta sequence allowed)')
    parser.add_argument('mutated', type=file, help='Fasta file of SNPs')

    args = parser.parse_args()
    plot_snps(**vars(args))
