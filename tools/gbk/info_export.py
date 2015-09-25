#!/usr/bin/env python
import argparse
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

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
    if isinstance(comment, list):
        return '%s\t%s' % (genome_id, '\t'.join(comment))
    return '%s\t%s' % (genome_id, comment)


def extract_metadata(genbank_file, section):
    with open(genbank_file, 'r') as handle:
        for record in SeqIO.parse(handle, "genbank"):
            output = ''
            if section in record.annotations:
                if section == 'references':
                    output += references_handler(record.annotations[section],
                                                 record.id)
                elif section in ['comment']:
                    output += comment_handler(record.annotations[section],
                                              record.id)
                elif section in ['source', 'taxonomy', 'date', 'organism']:
                    output += tabular_handler(record.annotations[section],
                                              record.id)
                else:
                    log.error("Cannot handle section type %s. "
                              "Please submit a feature request", section)
            yield output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Genbank Metadata Export')
    parser.add_argument('genbank_files', nargs='+', type=file, help='Genbank file')
    parser.add_argument('--section', type=str, help='Section to export', default='taxonomy')

    args = parser.parse_args()

    for gbk_file in args.genbank_files:
        for metadata in extract_metadata(gbk_file, args.section):
            print metadata
