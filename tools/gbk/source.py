#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
import logging
import cpt
from Bio import SeqIO
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def extract_host_info(host_string):
    source_host = host_string.split(' ')
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
        species = ""

    if len(source_host) > 0:
        genus = source_host[0]
    else:
        genus = ""
    return (id, genus, species, strain, other)


def phage_source(genbank_files=None, **kwargs):
    for genbank_file in genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            id = record.id

            source_feats = [x for x in record.features if x.type == 'source']

            # Provide a default value from a parsing attempt at the name
            (host, phage) = cpt.phage_name_parser(record.description)
            if host is not None:
                ret = (id, host)
            else:
                ret = (id, '')

            # If there isn't a source feature, return that immediately
            if len(source_feats) == 0:
                yield ret
            elif len(source_feats) == 1:
                # Otherwise look through any and all source features for a /host qualifier
                if 'host' in source_feats[0].qualifiers:
                    ret = extract_host_info(source_feats[0].qualifiers['host'][0])
                yield ret
            else:
                log.warning("Record with multiple source features. Unlikely to "
                            "be well behaved. Please email the developer with "
                            "this record %s", record.id)
                if 'host' in source_feats[0].qualifiers:
                    ret = extract_host_info(source_feats[0].qualifiers['host'][0])
                yield ret

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract lineage information from genbank files')
    parser.add_argument('genbank_files', type=file, nargs='+', help='Genbank file')
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()
    print '\t'.join(['# ID', 'Host Genus', 'Host Species', 'Host Strain', 'Other'])

    for line in phage_source(**vars(args)):
        print '\t'.join(map(str, line))
