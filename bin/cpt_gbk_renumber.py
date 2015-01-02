#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
Gene Renumbering Tool
=====================

Renumber genes in a genome
"""


def renumber_genes(gbk_file=None, tag_to_update="locus_tag",
                   string_prefix="display_id", leading_zeros=3):

    from Bio import SeqIO
    records = list(SeqIO.parse(gbk_file, "genbank"))
    delta = {
        'Sheet1': {
            'header': ['old', 'new'],
            'data': []
        }
    }
    for i in range(len(records)):
        if string_prefix == 'display_id':
            format_string = records[i].id + '_%0' + str(leading_zeros) + 'd'
        else:
            format_string = string_prefix + '%0' + str(leading_zeros) + 'd'

        f_cds = [f for f in records[i].features if f.type == 'CDS']
        f_rbs = [f for f in records[i].features if f.type == 'RBS']
        f_gene = [f for f in records[i].features if f.type == 'gene']
        f_oth = [f for f in records[i].features if f.type not in ['CDS', 'RBS',
                                                                  'gene']]

        # Make sure we've hit every RBS and gene
        h_rbs = {}
        h_gene = {}
        for cds in f_cds:
            # If there's an associated gene feature, it will share a stop codon
            associated_genes = [f for f in f_gene if f.location.end ==
                                cds.location.end]
            # If there's an RBS it'll be upstream a bit.
            if cds.location.strand > 0:
                associated_rbss = [f for f in f_rbs if f.location.end <
                                   cds.location.start and f.location.end >
                                   cds.location.start - 30]
            else:
                associated_rbss = [f for f in f_rbs if f.location.start >
                                   cds.location.end and f.location.start <
                                   cds.location.end + 30]
            tmp_result = [cds]
            if len(associated_genes) > 0:
                h_gene[f_gene.index(associated_genes[0])] = True
                tmp_result.append(associated_genes[0])
            if len(associated_rbss) > 0:
                h_rbs[f_rbs.index(associated_rbss[0])] = True
                tmp_result.append(associated_rbss[0])
            # We choose to append to f_other as that has all features not
            # already accessed. It may mean that some gene/RBS features are
            # missed if they aren't detected here, which we'll need to handle.
            f_oth.append(tmp_result)

        def delta_old(feature):
            # First part of delta entry, old name
            val = '%s %s %s' % (f.location.start, f.location.end,
                                f.location.strand,)
            if tag_to_update in f.qualifiers:
                val = f.qualifiers[tag_to_update][0]
            return val

        # Why must these people start at 1
        index = 1
        # Because we have to modify the array we work on a separate one
        clean_features = []
        for feature in f_oth:
            # Lists are groups of cds+gene?+rbs? that need to be renumbered.
            if isinstance(feature, list):
                 #Renumber all features in list
                for f in feature:
                    delta_entry = [delta_old(f)]
                    f.qualifiers[tag_to_update] = [format_string % index]
                    clean_features.append(f)
                    delta_entry.append(format_string % index)
                    delta['Sheet1']['data'].append(delta_entry)
            else:
                delta_entry = [delta_old(f)]
                print f
                feature.qualifiers[tag_to_update] = [format_string % index]
                clean_features.append(feature)
                delta_entry.append(format_string % index)
                delta['Sheet1']['data'].append(delta_entry)
            index += 1

        # Update all features
        records[i].features = clean_features
    return (records, delta)


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['file', 'Genbank file renumber genes on',
             {'required': True, 'validate': 'File/Input'}],
            ['tag_to_update', 'Which tag is used to store gene numbers',
             {'validate': 'String', 'default': 'locus_tag'}],
            ['string_prefix', 'A string to use as a prefix for the numbering. Will be used as XXXXXXNNN where XXXXXX is the string and NNN is a numerical identifier. Using "display_id" has special meaning, it will use the genome\'s name/accession number',
             {'validate': 'String', 'default': 'display_id'}],
            ['leading_zeros', 'Number of leading zeros/padding',
             {'validate': 'Int', 'default': 3}]
        ],
        outputs=[
            [
                'genbank',
                'Renumbered Genbank File',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'renumbered',
                    'data_format': 'genomic/annotated',
                    'default_format': 'Genbank',
                }
            ],
            [
                'change_table',
                'Table of updated gene names',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'renumbered',
                    'data_format': 'text/tabular',
                    'default_format': 'CSV',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.genbank.RelabelTags',
            'appname': 'Relabel Genbank Genes',
            'appvers': '1.95',
            'appdesc': 're-labels genbank tags according to rules',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    (result, delta) = renumber_genes(gbk_file=options['file'],
                                     tag_to_update=options['tag_to_update'],
                                     string_prefix=options['string_prefix'],
                                     leading_zeros=options['leading_zeros'])
    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='genbank', GGO=opts)
    of.CRR(data=result)
    of = OutputFiles(name='change_table', GGO=opts)
    of.CRR(data=delta)
