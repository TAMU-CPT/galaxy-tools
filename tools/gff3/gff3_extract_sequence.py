#!/usr/bin/env python
import sys
import argparse
import logging
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gff3 import feature_lambda, feature_test_type, get_id
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def main(fasta, gff3, feature_filter=None, nodesc=False):

    if feature_filter == 'nice_cds':
        from gff2gb import gff3_to_genbank
        for rec in gff3_to_genbank(gff3, fasta):
            for feat in rec.features:
                if feat.type != 'CDS':
                    continue

                if nodesc:
                    description = ''
                else:
                    feat.qualifiers['ID'] = [feat._ID]
                    product = feat.qualifiers.get('product', '')
                    description = '{1} [Location={0.location};ID={0.qualifiers[ID][0]}]'.format(feat, product)

                yield [
                    SeqRecord(
                        feat.extract(rec).seq,
                        id=feat.qualifiers.get('locus_tag', get_id(feat)),
                        description=description
                    )
                ]

    else:
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        for rec in GFF.parse(gff3, base_dict=seq_dict):
            for feat in feature_lambda(
                rec.features,
                feature_test_type,
                {'type': feature_filter},
                subfeatures=False
            ):
                id = feat.id
                if len(id) == 0:
                    id = get_id(feat)

                if nodesc:
                    description = ''
                else:
                    important_data = {
                        'Location': feat.location,
                    }
                    if 'Name' in feat.qualifiers:
                        important_data['Name'] = feat.qualifiers.get('Name', [''])[0]

                    description = '[{}]'.format(
                        ';'.join([
                            '{key}={value}'.format(key=k, value=v) for (k, v) in important_data.iteritems()
                        ])
                    )

                yield [
                    SeqRecord(
                        feat.extract(rec).seq,
                        id=id,
                        description=description
                    )
                ]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export corresponding sequence in genome from GFF3', epilog="")
    parser.add_argument('fasta', type=file, help='Fasta Genome')
    parser.add_argument('gff3', help='GFF3 File')
    parser.add_argument('--feature_filter', default=None, help='Filter for specific feature types')
    parser.add_argument('--nodesc', action='store_true', help='Strip description field off')
    args = parser.parse_args()

    for seq in main(**vars(args)):
        SeqIO.write(seq, sys.stdout, 'fasta')
