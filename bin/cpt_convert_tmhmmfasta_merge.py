#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import sys
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
Merge TMHMM and Fasta file to Tabular Output
============================================

"""

def parse_tmhmm(key, tmd_list):
    n = tmd_list[0][2]
    c = tmd_list[-1][2]
    tmds = 0
    for line_item in tmd_list:
        if line_item[2] == 'TMhelix':
            tmds += 1
    return [key, n, c, tmds]


def extract_features(fa=None, tmhmm=None):
    output = {
        'Sheet1': {
            'header': ['Fasta ID', 'N', 'C', 'TMDs'],
            'data': [],
        }
    }
    from Bio import SeqIO
    records = list(SeqIO.parse(fa, "fasta"))
    id_list = []
    for i in range(len(records)):
        id_list.append(records[i].id)

    tmhmm_data = {}
    for line in (line for line in tmhmm.readlines() if not
                 line.startswith('#')):
        parts = line.split('\t')
        if parts[0] not in tmhmm_data:
            tmhmm_data[parts[0]] = []
        tmhmm_data[parts[0]].append(parts)

    for key in id_list:
        if key in tmhmm_data:
            output['Sheet1']['data'].append(parse_tmhmm(key, tmhmm_data[key]))
        else:
            output['Sheet1']['data'].append([key, 'n/a', 'n/a', '0'])

    return output


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['fasta', 'Fasta File', {'required': True, 'validate':
                                                'File/Input'}],
            ['tmhmm', 'TMHMM Output', {'required': True, 'validate':
                                                'File/Input'}],
        ],
        outputs=[
            [
                'results',
                'Tabular Results',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'export',
                    'data_format': 'text/tabular',
                    'default_format': 'TSV_U',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.generic.TMHMM_merger',
            'appname': 'Reformat TMHMM results',
            'appvers': '1.94',
            'appdesc': 'convert to easily parseable tabular format',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    result = extract_features(fa=options['fasta'],tmhmm=options['tmhmm'])
    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='results', GGO=opts)
    of.CRR(data=result)
