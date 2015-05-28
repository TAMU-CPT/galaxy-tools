#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import sys
import logging
logging.basicConfig(level=logging.INFO)


def mga_to_gff3(mga_orf_table=None):
    if mga_orf_table is None:
        raise ValueError("Must specify orf file")

    orfs = ['##gff-version 3']
    for line in mga_orf_table.readlines():
        if not line.startswith('#'):
            (gene_id, start, end, strand, phase, complete, score, model,
                rbs_start, rbs_end, rbs_score) = line.strip().split('\t')

            if rbs_start != '-':
                rbs_line = [
                    'seq',
                    'MGA',
                    'Shine_Dalgarno_sequence',
                    rbs_start,
                    rbs_end,
                    '.',
                    strand,
                    '.',
                    'ID=rbs-%s' % gene_id,
                ]
                orfs.append('\t'.join([str(x) for x in rbs_line]))

            gff_line = [
                'seq',  # seqid
                'MGA',  # source
                'CDS',  # type
                start,  # start
                end,  # end
                score,  # score
                strand,  # strand
                '.',  # phase
                'ID=%s' % gene_id,  # attr
            ]
            orfs.append('\t'.join([str(x) for x in gff_line]))
    return "\n".join(orfs)


__doc__ = """
Convert MetaGeneAnnotator Table to GFF3
=======================================
"""

if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['file', 'MGA Output', {'required': True, 'validate':
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
            'appid': 'edu.tamu.cpt.generic.mga-to-gff3',
            'appname': 'MGA to GFF3',
            'appvers': '1.0',
            'appdesc': 'convert formats',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    result = mga_to_gff3(mga_orf_table=options['file'])

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='data', GGO=opts)
    of.CRR(data=result)
