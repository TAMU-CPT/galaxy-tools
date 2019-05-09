#!/usr/bin/env python
import argparse
from Bio import SeqIO
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name='charges')

HTML_HEADER = '<html><head><title>Charges Report</title></head><body>'
HTML_FOOTER = '</body></html>'


def charges_html(fasta, aa, fgColor, bgColor, width=120):
    colour_scheme = zip([x.upper() for x in aa], bgColor, fgColor)

    # CSS and header styling
    css = """<style type="text/css">
    .list li { list-style: none; margin:10px }
    .info { float:left; width:20px }
    pre { font-size:1.3em }
    """
    info = '<h1>Charges</h1><h3>Legend</h3><ul class="list">'
    for group in colour_scheme:
        css += '.%s{ background: %s; color: %s}\n' % group
        info += '<li><span class="%s" style="padding:5px">%s</span></li>\n' % (group[0], group[0])
    css += '</style>'
    info += '</ul>'

    # Pre-calculate, so we can use for testing 'in'
    match_list = [group[0] for group in colour_scheme]

    page = ''
    # Parse sequences from fasta file
    for record in SeqIO.parse(fasta, "fasta"):
        page += '<pre><h3>&gt;%s %s</h3>\n' % (record.id, record.description)
        seq = list(str(record.seq).upper())

        idx = 0
        for i in range(0, len(seq), width):
            line_charges = []
            line_residues = seq[i:i + width]
            line_numbers = []

            for char in range(len(line_residues)):
                if line_residues[char] in 'KRkr':
                    line_charges.append('+')
                elif line_residues[char] in 'DEde':
                    line_charges.append('-')
                else:
                    line_charges.append(' ')

                # Could be swapped out for math with i+char...
                idx += 1
                if idx % 10 == 0:
                    line_numbers.append('%10s' % idx)

                # Replace with <span>
                for m in match_list:
                    if line_residues[char].upper() in m:
                        line_residues[char] = '<span class="%s">%s</span>' % (m, line_residues[char])

            page += ''.join(line_charges) + '\n'
            page += ''.join(line_residues) + '\n'
            page += ''.join(line_numbers) + '\n'
            page += '\n'
        page += '</pre>'
    return HTML_HEADER + css + info + page + HTML_FOOTER


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Top related genomes')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='Fasta protein file')
    parser.add_argument('--width', type=int, help='Plot width', default=120)
    parser.add_argument('--aa', nargs='+')
    parser.add_argument('--fgColor', nargs='+')
    parser.add_argument('--bgColor', nargs='+')
    args = parser.parse_args()

    print(charges_html(**vars(args)))
