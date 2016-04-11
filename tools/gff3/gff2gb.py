"""Convert a GFF and associated FASTA file into GenBank format.

Usage:
gff_to_genbank.py <GFF annotation file> <FASTA sequence file>
"""
import argparse
import sys
import re

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from BCBio import GFF
from gff3 import feature_lambda
default_name = re.compile("^gene_(\d+)$")


def rename_key(ds, k_f, k_t):
    """Rename a key in a dictionary and return it, FP style"""
    ds[k_t] = ds[k_f]
    del ds[k_f]
    return ds


def gff3_to_genbank(gff_file, fasta_file):
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)

    def test_true(feature, **kwargs):
        return True

    for record in gff_iter:
        full_feats = []
        features = record.features
        for feature in features:
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

        # We'll work on a separate copy of features to avoid modifying a list
        # we're iterating over
        replacement_feats = []
        # We'll re-do gene numbering while we're at it
        fid = 0
        # Renumbering requires sorting
        for feat in sorted(
                feature_lambda(features, test_true, {}, subfeatures=True),
                # Based on feature location
                key=lambda x:int(x.location.start)):

            # Our modifications only involve genes
            if feat.type == 'gene':
                fid += 1
                # Which have mRNAs we'll drop later
                for mRNA in feat.sub_features:
                    # And some exons below that
                    for sf in mRNA.sub_features:
                        # We set a locus_tag
                        sf.qualifiers['locus_tag'] = 'gene_%03d' % fid
                        # and if it is the CDS (not the RBS)
                        if sf.type == 'exon' and len(sf) > 15:
                            # We copy over the parent gene's data. Argh, apollo.
                            sf.qualifiers.update(feat.qualifiers)
                # Wipe out the parent gene's data, leaving only a locus_tag
                feat.qualifiers = {
                    'locus_tag': 'gene_%03d' % fid,
                }
            elif feat.type in ('mRNA', 'exon'):
                # already been handled. EWW.
                continue

            # Patch our features back in (even if they're non-gene features)
            replacement_feats.append(feat)
        features = replacement_feats

        # Meat of our modifications
        for flat_feat in feature_lambda(
                features, test_true, {}, subfeatures=True):

            # We use the full GO term, but it should be less than that.
            if flat_feat.type == 'Shine_Dalgarno_sequence':
                flat_feat.type = 'RBS'

            # Drop mRNAs, apollo crap, useless CDSs
            if flat_feat.type in ('mRNA', 'non_canonical_three_prime_splice_site', 'non_canonical_five_prime_splice_site', 'CDS'):
                continue

            # Try and figure out a name. We gave conflicting instructions, so
            # this isn't as trivial as it should be.
            protein_product = flat_feat.qualifiers.get('product', flat_feat.qualifiers.get('Product', [None]))[0]
            if 'Name' not in flat_feat.qualifiers:
                protein_name = '__none__'
            else:
                protein_name = flat_feat.qualifiers['Name'][0]
            if protein_product is None:
                protein_product = protein_name

            # If the feature is an exon (there will be two (+?))
            if flat_feat.type == 'exon':
                # If the exon is SHORT that means it's the
                # Shine_Dalgarno_sequence, otherwise it's the CDS
                if len(flat_feat) <= 15:
                    flat_feat.type = 'RBS'
                else:
                    flat_feat.type = 'CDS'

                if 'Name' in flat_feat.qualifiers:
                    del flat_feat.qualifiers['Name']
            # if 'ID' in flat_feat.qualifiers:
                # flat_feat.qualifiers['locus_tag'] = flat_feat.qualifiers['ID']

            for x in ('source', 'phase', 'Parent', 'ID', 'owner',
                    'date_creation', 'date_last_modified'):
                if x in flat_feat.qualifiers:
                    del flat_feat.qualifiers[x]

            # Add product tag
            if flat_feat.type == 'CDS':
                flat_feat.qualifiers['Product'] = protein_product
                # Wipes out any 'product' key
                flat_feat.qualifiers = rename_key(flat_feat.qualifiers, 'Product', 'product')

            # In genbank format, note is lower case.
            flat_feat.qualifiers = rename_key(flat_feat.qualifiers, 'Note', 'note')

            # Append the feature
            full_feats.append(flat_feat)

        # Update our features
        record.features = full_feats
        yield record

if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Convert gff3 to gbk')
    parser.add_argument('gff_file', type=file, help='GFF3 file')
    parser.add_argument('fasta_file', type=file, help='Fasta Input')
    args = parser.parse_args()

    for record in gff3_to_genbank(**vars(args)):
        SeqIO.write([record], sys.stdout, "genbank")
