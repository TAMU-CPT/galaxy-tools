#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import numpy
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
Find top related genomes
========================

Phage genomes share a huge percentage similarity, this tool makes the process
for figuring out top related genomes signifcantly easier
"""


def top_related(parent=None, blast_results=None, method='n'):
    hits = {}
    for line in blast_results.readlines():
        hit = line.split('\t')[1].split('|')[0]
        if hit in hits:
            hits[hit] += 1
        else:
            hits[hit] = 1
    # Top results
    top_accessions = []
    # Allow different result process methods
    if method == 'sd':
        std = numpy.std(hits.values())
        mean = numpy.mean(hits.values())
        for key, value in hits.iteritems():
            if value > (mean + 2 * std):
                top_accessions.append(key)
    else:
        for key, value in reversed(sorted(hits.iteritems(), key=lambda (k, v):
                                          (v, k))):
            if len(top_accessions) > 4:
                break
            else:
                top_accessions.append(key)

    print top_accessions

    from Bio import SeqIO
    records = list(SeqIO.parse(parent, "genbank"))
    top_records = []
    for i in range(len(records)):
        if records[i].id in top_accessions:
            top_records.append(records[i])
    # Return top
    return (top_records, top_accessions)


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['db_gbk', 'Original Genbank file used in Blast database',
             {'required': True, 'validate': 'File/Input'}],
            ['blast', 'Blast results', {'required': True, 'validate':
                                        'File/Input'}],
            ['method', 'Method to select "top" results',
             {'required': True, 'validate': 'Option', 'default': 'n',
              'options': {'n': 'Top 5 results', 'sd':
                          'Statistical significance'}}],
        ],
        outputs=[
            [
                'top_genomes',
                'Re-export of top matched genomes for use in downstream tools',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'top_genomes',
                    'data_format': 'genomic/annotated',
                    'default_format': 'Genbank',
                }
            ],
            [
                'accession_list',
                'Accession numbers of top matched genomes',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'top_accessions',
                    'data_format': 'text/plain',
                    'default_format': 'TXT',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.genbank.TopRelated',
            'appname': 'Top Related Genomes',
            'appvers': '1.94',
            'appdesc': 'finds top related sequences from blast results',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    (gbk, txt) = top_related(parent=options['db_gbk'],
                             blast_results=options['blast'],
                             method=options['method'])

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='top_genomes', GGO=opts)
    of.CRR(data=gbk)
    of2 = OutputFiles(name='accession_list', GGO=opts)
    of2.CRR(data='\n'.join(txt))
