#!/usr/bin/env python
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
from gff3 import feature_lambda, wa_unified_product_name, is_uuid, feature_test_type
default_name = re.compile("^gene_(\d+)$")


def rename_key(ds, k_f, k_t):
    """Rename a key in a dictionary and return it, FP style"""
    # If they key is not in the dictionary, just return immediately
    if k_f not in ds:
        return ds

    # Otherwise, we check if the target key is in there
    if k_t in ds:
        # If it is, we need to append
        ds[k_t] += ds[k_f]
    else:
        # if not, we can just set.
        ds[k_t] = ds[k_f]

    # Remove source
    del ds[k_f]
    return ds


def gff3_to_genbank(gff_file, fasta_file):
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)

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

        # We'll work on a separate copy of features to avoid modifying a list
        # we're iterating over
        replacement_feats = []
        # We'll re-do gene numbering while we're at it
        fid = 0

        # These are NON-GENE features (maybe terminators? etc?)
        for feature in feature_lambda(record.features, feature_test_type, {'type': 'gene'}, subfeatures=False, invert=True, recurse=False):
            if feature.type == 'terminator':
                replacement_feats.append(feature)

        # Renumbering requires sorting
        for feature in sorted(feature_lambda(record.features, feature_test_type, {'type': 'gene'}, subfeatures=True),
                              key=lambda x: int(x.location.start)):
            # Our modifications only involve genes
            fid += 1
            # Which have mRNAs we'll drop later
            for mRNA in feature.sub_features:
                # And some exons below that
                sf_replacement = []
                for sf in mRNA.sub_features:
                    # We set a locus_tag on ALL sub features
                    sf.qualifiers['locus_tag'] = 'gene_%03d' % fid
                    # Remove Names which are UUIDs
                    if is_uuid(sf.qualifiers['Name'][0]):
                        del sf.qualifiers['Name']

                    # If it is the RBS exon (mis-labelled by apollo as 'exon')
                    if sf.type == 'exon' and len(sf) < 10:
                        sf.type = 'Shine_Dalgarno_sequence'
                        sf_replacement.append(sf)
                    # and if it is the CDS
                    elif sf.type == 'CDS':
                        # Update CDS qualifiers with all info that was on parent
                        sf.qualifiers.update(feature.qualifiers)
                        sf_replacement.append(sf)
                    elif sf.type == 'tRNA':
                        sf.qualifiers.update(feature.qualifiers)
                        sf_replacement.append(sf)

                # Replace the subfeatures on the mRNA
                mRNA.sub_features = sf_replacement
            # Wipe out the parent gene's data, leaving only a locus_tag
            feature.qualifiers = {
                'locus_tag': 'gene_%03d' % fid,
            }

            # Patch our features back in (even if they're non-gene features)
            replacement_feats.append(feature)

        # Meat of our modifications
        for flat_feat in feature_lambda(
                replacement_feats, lambda x: True, {}, subfeatures=True):

            # We use the full GO term, but it should be less than that.
            if flat_feat.type == 'Shine_Dalgarno_sequence':
                flat_feat.type = 'RBS'

            # Drop mRNAs, apollo crap, useless CDSs
            if flat_feat.type in ('mRNA', 'non_canonical_three_prime_splice_site', 'non_canonical_five_prime_splice_site'):
                continue

            if flat_feat.type == 'CDS' and len(flat_feat) < 10:
                # Another RBS mistake
                continue

            # Try and figure out a name. We gave conflicting instructions, so
            # this isn't as trivial as it should be.
            protein_product = wa_unified_product_name(flat_feat)

            for x in ('source', 'phase', 'Parent', 'ID', 'owner',
                      'date_creation', 'date_last_modified'):
                if x in flat_feat.qualifiers:
                    if x == 'ID':
                        flat_feat._ID = flat_feat.qualifiers['ID']
                    del flat_feat.qualifiers[x]

            # Add product tag
            if flat_feat.type == 'CDS':
                flat_feat.qualifiers['Product'] = protein_product
                # Wipes out any 'product' key
                flat_feat.qualifiers = rename_key(flat_feat.qualifiers, 'Product', 'product')

            elif flat_feat.type == 'terminator':
                flat_feat.type = 'regulatory'
                flat_feat.qualifiers = {
                    'regulatory_class': 'terminator',
                }

            # In genbank format, note is lower case.
            flat_feat.qualifiers = rename_key(flat_feat.qualifiers, 'Note', 'note')
            flat_feat.qualifiers = rename_key(flat_feat.qualifiers, 'description', 'note')
            flat_feat.qualifiers = rename_key(flat_feat.qualifiers, 'protein', 'note')
            flat_feat.qualifiers = rename_key(flat_feat.qualifiers, 'Dbxref', 'db_xref')
            if 'Name' in flat_feat.qualifiers:
                del flat_feat.qualifiers['Name']

            # Append the feature
            full_feats.append(flat_feat)

        # Update our features
        record.features = full_feats
        # Strip off record names that would cause crashes.
        record.name = record.name[0:16]
        yield record


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Convert gff3 to gbk')
    parser.add_argument('gff_file', type=argparse.FileType("r"), help='GFF3 file')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), help='Fasta Input')
    args = parser.parse_args()

    for record in gff3_to_genbank(**vars(args)):
        SeqIO.write([record], sys.stdout, "genbank")
