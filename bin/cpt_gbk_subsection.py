#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import logging
logging.basicConfig(level=logging.INFO)
import copy


__doc__ = """
Gene Renumbering Tool
=====================

Renumber genes in a genome
"""


def genbank_subsection(file=None, start=None, end=None, *args, **kwargs):
    from Bio import SeqIO
    records = list(SeqIO.parse(file, "genbank"))
    fixed_records = []
    for i in range(len(records)):
        fixed_records.append(records[i][start-1:end])
    return fixed_records

if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['file', 'Genbank file renumber genes on',
             {'required': True, 'validate': 'File/Input'}],
            ['start', 'Start of subsection',
             {'validate': 'Int', 'default': '1'}],
            ['end', 'End of subsection',
             {'validate': 'Int', 'default': '100'}],
        ],
        outputs=[
            [
                'subsection',
                'Subsection of Genbank File',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'subsection',
                    'data_format': 'genomic/annotated',
                    'default_format': 'Genbank',
                }
            ],
        ],
        defaults={
            'appid': 'edu.tamu.cpt.genbank.subsection',
            'appname': 'Subsection of Genbank File',
            'appvers': '0.1',
            'appdesc': 'take a Genbank file and produce a subsection retaining features',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    result = genbank_subsection(**options)

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='subsection', GGO=opts)
    of.CRR(data=result)
