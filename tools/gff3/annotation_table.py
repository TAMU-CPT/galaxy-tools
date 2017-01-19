#!/usr/bin/env python
import argparse
from gff2gb import gff3_to_genbank

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def to_annotation_table(seqrecord):
    # uid, length, strand, product, name, notes
    for feature in seqrecord.features:
        if feature.type not in ('CDS', 'tRNA', 'terminator'):
            continue

        print '\t'.join(map(str, (
            feature.id,
            feature.type,
            len(feature),
            '-' if feature.location.strand < 0 else '+',
            feature.location.start,
            feature.location.end,
            feature.qualifiers.get('product', feature.qualifiers['locus_tag']),
            feature.qualifiers['locus_tag'],
            '. '.join(feature.qualifiers.get('note', []))
        )))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='produce annotation table like product')
    parser.add_argument('gff_file', type=file, help='GFF3 file')
    parser.add_argument('fasta_file', type=file, help='Fasta Input')
    args = parser.parse_args()

    print '# ' + 'ID type length strand start end product locus_tag Notes'.replace(' ', '\t')
    for record in gff3_to_genbank(**vars(args)):
        print to_annotation_table(record)
