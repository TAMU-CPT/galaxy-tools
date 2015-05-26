#!/usr/bin/env python
import argparse
from Bio import SeqIO
from cpt_ncbi_efetch_gbk import name_parser


def phage_source(genbank_files=None, **kwargs):
    for genbank_file in genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            id = record.id
            if '.' in id:
                id = id.split('.')[0]

            source_feats = [x for x in record.features if x.type == 'source']

            # Provide a default value from a parsing attempt at the name
            (host, phage) = name_parser(record.description)
            if host is not None:
                ret = (id, host)
            else:
                ret = (id, '')

            # If there isn't a source feature, return that immediately
            if len(source_feats) == 0:
                yield ret
            else:
                # Otherwise look through any and all source features for a /host qualifier
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
                        ret = (id, genus, species, strain, other)
                yield ret

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract lineage information from genbank files')
    parser.add_argument('genbank_files', type=file, nargs='+', help='Genbank file')
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()
    print '\t'.join(['# ID', 'Host Genus', 'Host Species', 'Host Strain', 'Other'])

    for line in phage_source(**vars(args)):
        print '\t'.join(map(str, line))
