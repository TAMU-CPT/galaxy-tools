#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
Fasta Reopen
============

Reopen a fasta contig

For the example contig:

::

    >test
    aaaaaccccc

Reopening with 5, implying that we should reopen after the fifth base, we see:

::

    >test
    cccccaaaaa
"""


def reopen_contig(fasta_file=None, after=1000, **kwd):
    from Bio import SeqIO
    records = list(SeqIO.parse(fasta_file, "fasta"))
    output = []
    for i in range(len(records)):
        seq = records[i].seq
        records[i].seq = seq[after:] + seq[0:after]
        output.append(records[i])

    return output


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['fasta_file', 'Fasta file to reopen', {'required': True, 'validate':
                                              'File/Input'}],
            ['after', 'Reopen contig after this base (1-indexed)', {'required': True, 'validate': 'Int', 'default': '1000'}],
        ],
        outputs=[
            [
                'fasta',
                'Reopened fasta file',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'export',
                    'data_format': 'genomic/raw',
                    'default_format': 'Fasta',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.fasta.reopen',
            'appname': 'Reopen Fasta Sequences',
            'appvers': '1.0',
            'appdesc': 'at specific location',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    result = reopen_contig(**options)

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='fasta', GGO=opts)
    of.CRR(data=result)
