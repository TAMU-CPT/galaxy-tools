#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
Simple Primer Tool
==================

For the example contig:

::

    >test
    aaaaaNNNNNNNNNNNNccccc

Requesting a primer with a 5 base end overlap, you would get

::

    >test
    ccccc[ccc]aaaaa
"""


def generate_primers(fasta_file=None, end_overlap=500, **kwd):
    from Bio import SeqIO
    records = list(SeqIO.parse(fasta_file, "fasta"))
    output = ""
    for i in range(len(records)):
        seq = records[i].seq
        output += ">%s %s\n%s[%s]%s\n" % (
            records[i].id,
            records[i].description,
            seq[-end_overlap:],
            seq[-3:],
            seq[0:end_overlap]
        )

    return output


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['fasta_file', 'Fasta file to reopen', {'required': True, 'validate':
                                              'File/Input'}],
            ['end_overlap', 'End overlap', {'required': True, 'validate': 'Int', 'default': '500'}],
        ],
        outputs=[
            [
                'fasta',
                'Primer output',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'primers',
                    'data_format': 'text/plain',
                    'default_format': 'TXT',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.fasta.simpleprimer',
            'appname': 'Simple Primer Generator',
            'appvers': '1.0',
            'appdesc': '',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    result = generate_primers(**options)

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='fasta', GGO=opts)
    of.CRR(data=result)
