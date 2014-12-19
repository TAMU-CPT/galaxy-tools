#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import sys
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
GenBank Metadata Export
=======================

Export metadata from genbank files like taxonomy, citations, etc.
"""
bib_template = """
%% genome/ref: %s/%s
@Article{%s,
    Author="%s",
    Title="%s",
    Journal="%s",
}
"""


def references_handler(ref_list, genome_id):
    bib_list = ''

    i = 0
    for ref in ref_list:
        i += 1
        pmid = ref.pubmed_id
        if pmid == '':
            pmid = '%s.%s' % (genome_id, i)
        bib_list += bib_template % (genome_id, i, pmid, ref.authors, ref.title,
                                    ref.journal)
    return bib_list


def comment_handler(comment, genome_id):
    return '=== %s ===\n\n%s\n\n' % (genome_id, comment)


def tabular_handler(comment, genome_id):
    return '%s\t%s\n' % (genome_id, comment)


def extract_metadata(gbk_file=None, section=None):
    from Bio import SeqIO
    records = list(SeqIO.parse(gbk_file, "genbank"))
    output = ''
    for i in range(len(records)):
        if section in records[i].annotations:
            if section == 'references':
                output += references_handler(records[i].annotations[section],
                                             records[i].id)
            elif section in ['comment']:
                output += comment_handler(records[i].annotations[section],
                                          records[i].id)
            elif section in ['source', 'taxonomy', 'date', 'organism']:
                output += tabular_handler(records[i].annotations[section],
                                          records[i].id)
            else:
                #['comment', 'source', 'taxonomy', 'keywords', 'references',
                #'data_file_division', 'date', 'organism']
                print records[i].annotations.keys()
                output += "Don't know how to handle that section type"
    return output


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['file', 'Genbank file to filter', {'required': True, 'validate':
                                                'File/Input'}],
            ['section', 'Metadata elements to extract',
             {'required': True, 'validate': 'Option', 'multiple': False,
              'options': {'references': 'All references', 'taxonomy':
                          'Taxonomy', 'comment': "User entered genome comments (assembly data)",
                          'source': 'Genome source',
                          'keyworks': 'Genome keywords', 'date': 'Date',
                          'organism': 'Organism'}}],
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
            'appid': 'edu.tamu.cpt.genbank.MetadataExport',
            'appname': 'Genbank Metadata Export',
            'appvers': '1.94',
            'appdesc': 'exports metadata about genbank files (citations, taxonomy, etc)',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    result = extract_metadata(gbk_file=options['file'],
                              section=options['section'])
    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='data', GGO=opts)
    of.CRR(data=result)
