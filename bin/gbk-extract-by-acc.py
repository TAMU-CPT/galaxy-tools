#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
Extract Subset of Genbank File by Accessions
============================================

Sometimes you have a "master" database of genbank records, a single huge file
containing all of your genomes, then you do somet analysis which generates a
list of accession numbers, and you need to extract those from the main file for
downstream analysis.
"""


def extract_by_acc(parent=None, acc=None, strict=False):
    acc_list = []
    for line in acc.readlines():
        acc_list.append(line.strip())

    from Bio import SeqIO
    records = list(SeqIO.parse(parent, "genbank"))
    extracted = []
    for i in range(len(records)):
        acc = records[i].id
        acc = acc[0:acc.rindex('.')]
        if (strict and records[i].id in acc_list) or (not strict and acc in
           acc_list):
            extracted.append(records[i])
    # Return top
    return extracted


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['db_gbk', 'Genbank "database" to extract from',
             {'required': True, 'validate': 'File/Input'}],
            ['accession_list', 'List of Accessions to extract from' +
             'multi-record Genbank file', {'required': True, 'validate':
                                           'File/Input'}],
            ['strict', 'Accession numbers have versions, setting strict' +
             'indicates an exact version match is required'],
        ],
        outputs=[
            [
                'extracted',
                'Extracted genomes',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'extracted',
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
    gbk = extract_by_acc(parent=options['db_gbk'],
                         acc=options['accession_list'],
                         strict=options['strict'])

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='extracted', GGO=opts)
    of.CRR(data=gbk)
