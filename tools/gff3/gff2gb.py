"""Convert a GFF and associated FASTA file into GenBank format.

Usage:
gff_to_genbank.py <GFF annotation file> <FASTA sequence file>
"""
import argparse
import sys

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from BCBio import GFF
from gff3 import feature_lambda



def main(gff_file, fasta_file):
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)

    def test_true(feature, **kwargs):
        return True

    for record in gff_iter:
        full_feats = []
        for feature in record.features:
            if feature.type == 'region' and 'source' in feature.qualifiers and \
                    'GenBank' in feature.qualifiers['source']:
                feature.type = 'source'

                if 'comment1' in feature.qualifiers:
                    del feature.qualifiers['comment1']

                if 'Note' in feature.qualifiers:
                    record.annotations = feature.qualifiers
                    if len(feature.qualifiers['Note']) > 1:
                        record.annotations['comment'] = feature.qualifiers['Note'][1]
                    del feature.qualifiers['Note']

                if 'comment' in feature.qualifiers:
                    del feature.qualifiers['comment']

            for flat_feat in feature_lambda(
                    record.features, test_true, {}, subfeatures=True):
                if flat_feat.type == 'Shine_Dalgarno_sequence':
                    flat_feat.type = 'RBS'

                if 'ID' in flat_feat.qualifiers:
                    flat_feat.qualifiers['locus_tag'] = flat_feat.qualifiers['ID']

                for x in ('source', 'phase', 'Parent', 'ID'):
                    if x in flat_feat.qualifiers:
                        del flat_feat.qualifiers[x]
                full_feats.append(flat_feat)

        record.features = full_feats

        yield record

if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Convert gff3 to gbk')
    parser.add_argument('gff_file', type=file, help='GFF3 file')
    parser.add_argument('fasta_file', type=file, help='Fasta Input')
    args = parser.parse_args()

    for record in main(**vars(args)):
        SeqIO.write([record], sys.stdout, "genbank")
