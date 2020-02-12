##### Script will format a json that will house family members

import sys
import os
import json


endolysins = ['lysozyme',
              'endolysin',
              'lysin',
              'amidase',
              'muramidase',
              'muraminidase',
              'glycosidase',
              'endopeptidase'
              'muralytic',
              'glycosylase',
              'chitinase']

holins = ['holin',
          'antiholin',
          'hole',
          'pinholin',
          'lysis (protein)']


ECDs = ['Glyco_hydro_19',
        'Glyco_hydro_25',
        'Glyco_hydro_108',
        'Chitinase',
        'SLT',
        'Transglycosylase',
        'Glucosaminidase',
        'Amidase02_C',
        'Amidase_5',
        'Amidase_3',
        'Amidase_2',
        'NlpD',
        'VanY',
        'Peptidase_U40',
        'Peptidase_M15_3',
        'Peptidase_M15_4',
        'Peptidase_M23',
        'YkuD',
        'NLPC_P60',
        'Peptidase_C39_2',
        'CHAP',
        'DUF3597']

ECDs_accro = ['GH19',
            'GH25',
            'GH108',
            'GyH',
            'SLT',
            'TRANG',
            'GLUCO',
            'AMIO2-C',
            'AMI-5',
            'AMI-3',
            'AMI-2',
            'NLPD',
            'VANY',
            'PET-U40',
            'PET-15-3',
            'PET-15-4',
            'PET-M23',
            'YKUD',
            'NLPC-P60',
            'PET-C39-2',
            'CHAP',
            'DUF']

CBDs = ['PG_binding_3',
        'LysM',
        'SH3_3',
        'SH3_5',
        'PG_binding_1',
        'ChW',
        'ChW',
        'Cpl-7',
        'LGFP',
        'SH3-related',
        'FOG',
        'SH3b',
        'SPOR',
        'SLAP']

CBDs_accro = ['PG-3',
            'LYSM',
            'SH3-3',
            'SH3-5',
            'PG-1',
            'CHW',
            'CHW',
            'CPL7',
            'LGFP',
            'SH3-r',
            'FOG',
            'SH3b',
            'SPOR',
            'SLAP']

lysisCombo = {'endolysins':endolysins,'holins':holins,'ECDs':ECDs,'ECDs_accro':ECDs_accro,'CBDs':CBDs,'CBDs_accro':CBDs_accro}

filename = 'lysis-family.json'
with open('tools/proximity/data/'+filename, 'w') as j:
    json.dump(lysisCombo, j, indent='\t')