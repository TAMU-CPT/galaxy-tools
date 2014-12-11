#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import sys
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name='charges')

__doc__ = """
Charges
=======

"""

COLOUR_SCHEMES = {
    'default': {
        'polar': {
            # must be uppercase
            'match': 'HSQTNCY',
            'fg': 'black',
            'bg': 'white'
        },
        'hydrophobic': {
            'match': 'AVILMPFWG',
            'fg': 'black',
            'bg': '#999',
        },
        'charged': {
            'match': 'ERDK',
            'fg': 'white',
            'bg': 'black',
        }
    }
}
WIDTH = 120
HTML_HEADER = '<html><head><title>Charges Report</title></head><body>'
HTML_FOOTER = '</body></html>'

def charges_html(fasta=None, cs=None, **kwd):
    if cs not in COLOUR_SCHEMES:
        log.error("Colour scheme not known to script")
        sys.exit(1)

    # CSS and header styling
    css = """<style type="text/css">
    .list li { list-style: none }
    .info { float:left; width:20px }
    pre { font-size:1.3em }
    """
    info = '<ul class="list">'
    for group in COLOUR_SCHEMES[cs]:
        css += '.%(match)s{ background: %(bg)s; color: %(fg)s}\n' % COLOUR_SCHEMES[cs][group]
        info += '<li>&nbsp;-&nbsp;%(match)s<span class="info" style="background: %(bg)s">&nbsp;</span><span class="info" style="background: %(fg)s">&nbsp;</span></li>\n' % COLOUR_SCHEMES[cs][group]
    css += '</style>'
    info += '</ul>'

    # Pre-calculate, so we can use for testing 'in'
    match_list = [COLOUR_SCHEMES[cs][x]['match'] for x in COLOUR_SCHEMES[cs]]

    page = ''
    # Parse sequences from fasta file
    from Bio import SeqIO
    for record in SeqIO.parse(fasta, "fasta"):
        page += '<pre><h3>&gt;%s %s</h3>\n' % (record.id, record.description)
        seq = list(str(record.seq).upper())


        idx = 0
        for i in range(0, len(seq), WIDTH):
            line_charges = []
            line_residues = seq[i:i+WIDTH]
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
                    if line_residues[char] in m:
                        line_residues[char] = '<span class="%s">%s</span>' % (m, line_residues[char])


            page += ''.join(line_charges) + '\n'
            page += ''.join(line_residues) + '\n'
            page += ''.join(line_numbers) + '\n'
            page += '\n'
        page += '</pre>'
    return HTML_HEADER + css + info + page + HTML_FOOTER


def passthrough(cb):
    opts = GGO(
        options=[
            ['fasta', 'Fasta protein file',
             {'required': True, 'validate': 'File/Input'}],
            ['cs', 'Colour scheme',
             {'required': True, 'validate': 'Option', 'options':
              {
                  'default': 'Deafult hyrophobic/polar/charged'
              }
              }]
        ],
        outputs=[
            [
                'html',
                'HTML Report',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'charges',
                    'data_format': 'text/html',
                    'default_format': 'HTML',
                }
            ],
        ],
        defaults={
            'appid': 'edu.tamu.cpt.tools.charges',
            'appname': 'Charges',
            'appvers': '1.0',
            'appdesc': 'colour sequences based on rules',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    return (opts, options, cb(**options))


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    (opts, options, result) = passthrough(charges_html)

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='html', GGO=opts)
    of.CRR(data=result)
