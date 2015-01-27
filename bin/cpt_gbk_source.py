#!/usr/bin/env python
import argparse
import logging
from cpt_ncbi_efetch_gbk import name_parser
logging.basicConfig(level=logging.INFO)


__doc__ = """
Gene Renumbering Tool
=====================

Renumber genes in a genome
"""


def phage_source(genbank_files=None, **kwargs):
    from Bio import SeqIO
    rows = []
    for genbank_file in genbank_files:
        records = list(SeqIO.parse(genbank_file, "genbank"))
        for record in records:
            id = record.id
            if '.' in id:
                id = id.split('.')[0]

            source_feats = [x for x in record.features if x.type == 'source']
            if len(source_feats) == 0:
                (host, phage) = name_parser(record.description)
                rows.append([id, host])
            else:
                for source_feat in source_feats:
                    if 'host' in source_feat.qualifiers:
                        source_host = source_feat.qualifiers['host'][0].split(' ')
                        if len(source_host) > 3:
                            other = ' '.join(source_host[3:])
                        else:
                            other = ""

                        if len(source_host) > 2:
                            strain = source_host[2]
                        else:
                            strain = ""

                        if len(source_host) > 1:
                            species = source_host[1]
                        else:
                            species  = ""

                        if len(source_host) > 0:
                            genus = source_host[0]
                        else:
                            genus = ""

                        rows.append([id, genus, species, strain, other])
    return rows

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract lineage information from genbank files')
    parser.add_argument('genbank_files', type=file, nargs='+', help='Genbank file')
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()
    print '\t'.join(['# ID', 'Host Genus', 'Host Species', 'Host Strain', 'Other'])
    for line in phage_source(**vars(args)):
        print '\t'.join(line)
