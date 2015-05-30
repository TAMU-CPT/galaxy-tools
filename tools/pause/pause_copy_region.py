#!/usr/bin/env python
"""PAUSE: Plotter
"""
from Bio import SeqIO
import cpt_pause
from galaxygetopt.ggo import GalaxyGetOpt as GGO

def main(genomic_region=None, genomic_region_start=None, genomic_region_end=None, genome=None, **kwd):
    cut_start = 0
    cut_end = 0

    if genomic_region_start is not None and genomic_region_end is not None:
        cut_start = genomic_region_start
        cut_end = genomic_region_end
    elif genomic_region is not None:
        data = genomic_region.readlines()[1].split('\t')
        cut_start = data[0]
        cut_end = data[1]
    else:
        raise Exception("Must specify region")

    cut_start = int(cut_start)
    cut_end = int(cut_end)

    for record in SeqIO.parse(genome, "fasta") :
        reassemble = record.seq[cut_start:] + record.seq[0:cut_start] + record.seq[cut_start:cut_end]

        record.seq = reassemble
        return [record]

    return None

if __name__ == "__main__":
    opts = GGO(
        options=[
            ['genome', 'Input fasta genome',
             {'multiple': False, 'validate': 'File/Input', 'format': ['Fasta']}],
            ['genomic_region', 'Region of the genome which is a repeat',
             {'multiple': False, 'validate': 'File/Input'}],

            ['genomic_region_start', 'Start of repeat',
             {'multiple': False, 'validate': 'String'}],
            ['genomic_region_end', 'End of repeat',
             {'multiple': False, 'validate': 'String'}],
        ],
        outputs=[
            [
                'repeat_fa',
                'Fasta file with repeat added',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'repeat_fa',
                    'data_format': 'genomic/raw',
                    'default_format': 'Fasta',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.pause2.repeat_adder',
            'appname': 'PAUSE2 Repeat Adder',
            'appvers': '0.3',
            'appdesc': 'reopen+copy region based on a PAUSE result',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    output = main(**options)

    from galaxygetopt.outputfiles import OutputFiles
    off = OutputFiles(name='repeat_fa', GGO=opts)
    off.CRR(data=output)
