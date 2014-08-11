#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import sys
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
GenBank Feature Export
======================

Exports features from a GenBank file
"""


def get_id(feature=None, parent_prefix=None):
    result = ""
    if parent_prefix is not None:
        result += parent_prefix + '|'
    if 'locus_tag' in feature.qualifiers:
        result += feature.qualifiers['locus_tag'][0]
    elif 'gene' in feature.qualifiers:
        result += feature.qualifiers['gene'][0]
    elif 'product' in feature.qualifiers:
        result += feature.qualifiers['product'][0]
    else:
        result += '%s_%s_%s' % (feature.location.start, feature.location.end,
                                feature.location.strand)
    return result


def ensure_location_in_bounds(start=0, end=0, parent_length=0):
    # This prevents frameshift errors
    while start < 0:
        start += 3
    while end < 0:
        end += 3
    while start > parent_length:
        start -= 3
    while end > parent_length:
        end -= 3
    return (start, end)


def extract_features(gbk_file=None, tag='CDS', translate=False,
                     n_bases_upstream=0, n_bases_downstream=0,
                     strip_stops=False, tn_table=11):
    output = []
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    records = list(SeqIO.parse(gbk_file, "genbank"))
    for i in range(len(records)):
        for feature in records[i].features:
            if feature.type in tag:
                # Find new feature boundaries
                start = feature.location.start
                end = feature.location.end
                strand = feature.location.strand
                if n_bases_downstream == 0:
                    # Remove stops
                    if strip_stops:
                        if strand > 0:
                            end -= 3
                        else:
                            start += 3
                else:
                    # If we want extra on the end we cannot listen to
                    # stop_stripping requests
                    if strand > 0:
                        end += n_bases_downstream
                    else:
                        start -= n_bases_downstream
                # n_bases_upstream
                if strand > 0:
                    start -= n_bases_upstream
                else:
                    end += n_bases_upstream

                (start, end) = ensure_location_in_bounds(start=start, end=end,
                                                         parent_length=records[i].__len__)


                # Create our temp feature used to obtain correct portion of
                # genome
                tmp = SeqFeature(FeatureLocation(start, end, strand=strand),
                                 type='domain')
                # Translate the sequence
                if translate:
                    seq = tmp.extract(records[i].seq).translate(table=tn_table)
                else:
                    seq = tmp.extract(records[i].seq)
                output.append(SeqRecord(seq=seq, id=get_id(feature, parent_prefix=records[i].id),
                                        name=get_id(feature),
                                        description=get_id(feature)))
    return output


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['file', 'Genbank file to filter', {'required': True, 'validate':
                                                'File/Input'}],
            ['tag', 'Tag to extract', {'required': True, 'validate':
                                       'Genomic/Tag', 'multiple': True,
                                       'default': ['CDS']}],
            ['translate', 'Translate sequence during analysis'],
            ['translation_table_id', 'ID Number of tranlsation table to use',
             {'default': 11, 'validate': 'Int', 'required': True}],
            ['n_bases_upstream', 'Added N bases upstream to result',
             {'validate': 'Int', 'default': 0, 'min': 0}],
            ['n_bases_downstream', 'Added N bases downstream to result',
             {'validate': 'Int', 'default': 0, 'min': 0}],
            ['strip_stops', 'Remove all stop codons from translated proteins',
             {'validate': 'Flag'}],
        ],
        outputs=[
            [
                'fasta',
                'Fasta export of selected tag',
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
            'appid': 'edu.tamu.cpt.genbank.FeatureExport',
            'appname': 'Genbank Feature Export',
            'appvers': '1.94',
            'appdesc': 'Export features from a Genbank File',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    result = extract_features(gbk_file=options['file'], tag=options['tag'],
                              translate=options['translate'],
                              n_bases_upstream=options['n_bases_upstream'],
                              n_bases_downstream=options['n_bases_downstream'],
                              strip_stops=options['strip_stops'],
                              tn_table=options['translation_table_id'])
    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='fasta', GGO=opts)
    of.CRR(data=result)
