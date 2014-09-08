#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import sys
import logging
logging.basicConfig(level=logging.INFO)


def glimmer3_to_gff3(glimmer3_orf_table=None):
    if glimmer3_orf_table is None:
        raise ValueError("Must specify orf file")

    orfs = []
    for line in glimmer3_orf_table.readlines():
        if not line.startswith('>'):
            (g3i, g3s, g3e, g3p, g3score) = line.strip().split()
            strand = 1
            if g3s > g3e:
                strand = -1
                start = g3e
                end = g3s
            else:
                start = g3s
                end = g3e
            gff_line = [
                'glimmer',  # seqid
                'Glimmer3',  # source
                'CDS',  # type
                start,  # start
                end,  # end
                g3score,  # score
                strand,  # strand
                '.',  # phase
                '.',  # attr
            ]
            orfs.append('\t'.join([str(x) for x in gff_line]))
    return "\n".join(orfs)


__doc__ = """
Convert Glimmer3 Table to GFF3
==============================
"""

if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['file', 'Glimmer3 Output', {'required': True, 'validate':
                                         'File/Input'}],
        ],
        outputs=[
            [
                'data',
                'Exported data',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'export',
                    'data_format': 'text/plain',
                    'default_format': 'TXT',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.generic.glimmer3-to-gff3',
            'appname': 'Glimmer3 to GFF3',
            'appvers': '1.0',
            'appdesc': 'convert formats',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    result = glimmer3_to_gff3(glimmer3_orf_table=options['file'])

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='data', GGO=opts)
    of.CRR(data=result)
