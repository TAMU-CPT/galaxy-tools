#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
GenBank Feature Export
======================

Exports features from a GenBank file
"""


def extract_features(fasta_file=None, id="merged"):
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    sequence = ''
    ids = []
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for i in range(len(records)):
        ids.append(records[i].id)
        sequence += records[i].seq

    output = []
    output.append(SeqRecord(seq=sequence, id=id, description='Created from [%s...]' % (','.join(ids[0:10]))))
    return output


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['file', 'Fasta file to merge', {'required': True, 'validate':
                                             'File/Input'}],
            ['id', 'New fasta identifier for merged sequences', {'required': True, 'validate': 'String', 'default': 'merged'}],
        ],
        outputs=[
            [
                'fasta',
                'Concatenated fasta file',
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
            'appid': 'edu.tamu.cpt.genbank.fasta-merge',
            'appname': 'Concatenate Fasta Sequences',
            'appvers': '1.94',
            'appdesc': 'merge into single long sequence',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    result = extract_features(fasta_file=options['file'], id=options['id'])
    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='fasta', GGO=opts)
    of.CRR(data=result)
