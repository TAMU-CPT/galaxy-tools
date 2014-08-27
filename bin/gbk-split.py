#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
Split Multi-record Genbank Files
================================

Artemis will only open single-record genbank files.
"""


def split_gbks(parent=None):
    from Bio import SeqIO
    records = list(SeqIO.parse(parent, "genbank"))
    extracted = []
    for i in range(len(records)):
        extracted.append({
            'name': records[i].id,
            'gbk': records[i]
        })
    # Return top
    return extracted


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['gbk', 'Multi-record genbank file to split',
             {'required': True, 'validate': 'File/Input'}],
        ],
        outputs=[
            [
                'genbank',
                'Extracted genomes',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'split',
                    'data_format': 'genomic/annotated',
                    'default_format': 'Genbank',
                }
            ],
        ],
        defaults={
            'appid': 'edu.tamu.cpt.genbank.ExtractAccessions',
            'appname': 'Subset Genbank File: Accessions',
            'appvers': '1.94',
            'appdesc': 'given an accession list, extract a subset of the '
                       + 'genbank records into a new file.',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    gbks = split_gbks(parent=options['gbk'])

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='genbank', GGO=opts)
    for gbk in gbks:
        of.varCRR(data=gbk['gbk'], filename=gbk['name'])
