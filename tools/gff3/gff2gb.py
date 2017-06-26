#!/usr/bin/env python
"""Convert a GFF and associated FASTA file into GenBank format.

Usage:
gff_to_genbank.py <GFF annotation file> <FASTA sequence file>
"""
import argparse
import sys
import re
import copy
import itertools
import logging
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from BCBio import GFF
from gff3 import feature_lambda, wa_unified_product_name, is_uuid, \
    feature_test_type, fsort, feature_test_true, feature_test_quals
default_name = re.compile(r"^gene_(\d+)$")
logging.basicConfig(level=logging.INFO)


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
        yield handle_record(record)


def handle_non_gene_features(features):
    # These are NON-GENE features (maybe terminators? etc?)
    for feature in feature_lambda(features, feature_test_type, {'type': 'gene'}, subfeatures=False, invert=True, recurse=False):
        if feature.type in ('terminator', 'tRNA'):
            yield feature


def fminmax(feature):
    fmin = None
    fmax = None
    for sf in feature_lambda([feature], feature_test_true, {}, subfeatures=True):
        if fmin is None:
            fmin = sf.location.start
            fmax = sf.location.end
        if sf.location.start < fmin:
            fmin = sf.location.start
        if sf.location.end > fmax:
            fmax = sf.location.end
    return fmin, fmax


def fix_gene_boundaries(feature):
    # There is a frustrating bug in apollo whereby we have created gene
    # features which are LARGER than expected, but we cannot see this.
    # We only see a perfect sized gene + great SD together.
    #
    # So, we have this awful hack to clamp the location of the gene
    # feature to the contained mRNAs. This is good enough for now.
    fmin, fmax = fminmax(feature)
    if feature.location.strand > 0:
        feature.location = FeatureLocation(fmin, fmax, strand=1)
    else:
        feature.location = FeatureLocation(fmin, fmax, strand=-1)
    return feature


def fix_gene_qualifiers(name, feature, fid):
    for mRNA in feature.sub_features:
        mRNA.qualifiers['locus_tag'] = 'CPT_%s_%03d' % (name, fid)
        # And some exons below that
        sf_replacement = []
        for sf in mRNA.sub_features:
            # We set a locus_tag on ALL sub features
            sf.qualifiers['locus_tag'] = 'CPT_%s_%03d' % (name, fid)
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

        if mRNA.type == 'tRNA':
            mRNA.qualifiers['product'] = mRNA.qualifiers['Name']

        # Handle multiple child CDS features by merging them.
        # Replace the subfeatures on the mRNA
        mRNA.sub_features = merge_multi_cds(sf_replacement)
    return feature


def fix_frameshifted(features):
    logging.info("Fixing Frameshifted group: [%s]", str(features))
    genes = features
    # Find all mRNAs (plus reduce nested list into flattened one)
    mRNAs = sum([f.sub_features for f in genes], [])
    # Find all CDSs (plus reduce nested list into flattened one)
    cdss = sum([m.sub_features for m in mRNAs], [])
    # List to store the RBSs which we'll break apart + re-attach later.
    rbss = []
    # List to store all of the CDSs (i.e. cdss - rbss)
    cdss2 = []
    # Copy genes + clean out subfeatures. We'll re-use these constructs.
    fixed_features = copy.deepcopy(genes)
    for f in fixed_features:
        f.sub_features = []
    # Copy / empty out mRNAs
    fixed_mrnas = copy.deepcopy(mRNAs)
    for f in fixed_mrnas:
        f.sub_features = []
        f.qualifiers = {}
    # Fill rbss + cdss2
    for cds in cdss:
        if 'frameshift' in cds.qualifiers:
            del cds.qualifiers['frameshift']
        # Ignore short features, as those are RBSs
        if len(cds) < 15:
            rbss.append(cds)
            continue
        # Otherwise cdss.
        else:
            cdss2.append(cds)
    # Ok, now have cdss2 to deal with.
    other = []
    # Find the two with least value for distance between end / start (strand aware).
    # For every possible pair, we'll check their distance
    match_data = {}
    for (a, b) in itertools.permutations(cdss2, 2):
        if a.location.start < b.location.start:
            # A is downstream of B
            match_data[(a, b)] = b.location.start - a.location.end
        else:
            match_data[(a, b)] = a.location.start - b.location.end

    # Now we'll find the features which are closest in terms of start/end
    ((merge_a, merge_b), value) = max(match_data.items(), key=lambda kv: kv[1])
    # And get the non-matching features into other
    for f in cdss2:
        if f != merge_a and f != merge_b:
            other.append(f)
    # Back to the merge_a/b
    # With those, we'll merge them into one feature, and discard the other.
    merge_a.location = CompoundLocation([merge_a.location, merge_b.location])
    # The gene + RBSs should be identical and two/two.
    assert len(fixed_features) == 2
    # If not, we can just duplicate the RBS, doesn't matter.
    if len(rbss) != 2:
        rbss = [rbss[0], copy.deepcopy(rbss[0])]
    # Now re-construct.
    gene_0 = fixed_features[0]
    gene_1 = fixed_features[1]
    mRNA_0 = fixed_mrnas[0]
    mRNA_1 = fixed_mrnas[1]

    mRNA_0.sub_features = [rbss[0], merge_a]
    mRNA_1.sub_features = other + [rbss[1]]
    mRNA_0 = fix_gene_boundaries(mRNA_0)
    mRNA_1 = fix_gene_boundaries(mRNA_1)

    gene_0.sub_features = [mRNA_0]
    gene_1.sub_features = [mRNA_1]
    gene_0 = fix_gene_boundaries(gene_0)
    gene_1 = fix_gene_boundaries(gene_1)

    return fixed_features


def fix_frameshifts(features):
    # Collect all gene features where at least one subfeature has a
    # frameshift=??? annotation.
    def has_frameshift_qual(f):
        return len(list(feature_lambda(f.sub_features, feature_test_quals, {'frameshift': None}))) > 0

    def has_frameshift_qual_val(f, val):
        return len(list(feature_lambda(f.sub_features, feature_test_quals, {'frameshift': val}))) > 0

    def get_frameshift_qual(f):
        for f in feature_lambda(f.sub_features, feature_test_quals, {'frameshift': None}):
            return f.qualifiers['frameshift']

    to_frameshift = [x for x in features if x.type == 'gene' and has_frameshift_qual(x)]
    fixed = [x for x in features if x not in to_frameshift]

    frameshift_keys = set(sum(map(get_frameshift_qual, to_frameshift), []))
    for key in frameshift_keys:
        # Get features matching that key
        current = [x for x in to_frameshift if has_frameshift_qual_val(x, key)]
        # Fix them and append them
        fixed += fix_frameshifted(current)

    return fixed


def remove_useless_features(features):
    # Drop mRNAs, apollo crap, useless CDSs
    for f in features:
        if f.type in ('non_canonical_three_prime_splice_site',
                      'non_canonical_five_prime_splice_site',
                      'stop_codon_read_through', 'mRNA'):
            continue
        else:
            if f.type == 'CDS' and len(f) < 10:
                # Another RBS mistake
                continue
            # We use the full GO term, but it should be less than that.
            if f.type == 'Shine_Dalgarno_sequence':
                f.type = 'RBS'
            yield f


def merge_multi_cds(mRNA_sf):
    cdss = [x for x in mRNA_sf if x.type == 'CDS']
    non_cdss = [x for x in mRNA_sf if x.type != 'CDS']
    if len(cdss) <= 1:
        return non_cdss + cdss
    else:
        # Grab all locations, and sort them so we can work with them rationally.
        locations = sorted([x.location for x in cdss], key=lambda x: x.location.start)
        # Pick randomly a main CDS
        main_cds = cdss[0]
        # We'll merge the other CDSs into this one.
        main_cds.location = CompoundLocation(locations)
        return non_cdss + [main_cds]


def handle_record(record):
    full_feats = []
    for feature in fsort(record.features):
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
    replacement_feats += list(handle_non_gene_features(record.features))

    # Renumbering requires sorting
    fid = 0
    for feature in fsort(feature_lambda(record.features, feature_test_type, {'type': 'gene'}, subfeatures=True)):
        # Our modifications only involve genes
        fid += 1

        feature = fix_gene_boundaries(feature)
        # Which have mRNAs we'll drop later
        feature = fix_gene_qualifiers(record.id, feature, fid)

        # Wipe out the parent gene's data, leaving only a locus_tag
        feature.qualifiers = {
            'locus_tag': 'CPT_%s_%03d' % (record.id, fid),
        }

        # Patch our features back in (even if they're non-gene features)
        replacement_feats.append(feature)

    replacement_feats = fix_frameshifts(replacement_feats)
    flat_features = feature_lambda(replacement_feats, lambda x: True, {}, subfeatures=True)
    flat_features = remove_useless_features(flat_features)

    # Meat of our modifications
    for flat_feat in flat_features:

        # Try and figure out a name. We gave conflicting instructions, so
        # this isn't as trivial as it should be.
        protein_product = wa_unified_product_name(flat_feat)

        for x in ('source', 'phase', 'Parent', 'ID', 'owner', 'date_creation',
                  'date_last_modified'):
            if x in flat_feat.qualifiers:
                if x == 'ID':
                    flat_feat._ID = flat_feat.qualifiers['ID']
                del flat_feat.qualifiers[x]

        # Add product tag
        if flat_feat.type == 'CDS':
            flat_feat.qualifiers['product'] = [protein_product]
            if 'Product' in flat_feat.qualifiers:
                del flat_feat.qualifiers['Product']

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

        # more apollo nonsense
        if 'Manually set translation start' in flat_feat.qualifiers.get('note', []):
            flat_feat.qualifiers['note'].remove('Manually set translation start')

        # Append the feature
        full_feats.append(flat_feat)

    # Update our features
    record.features = fsort(full_feats)
    # Strip off record names that would cause crashes.
    record.name = record.name[0:16]
    return record


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Convert gff3 to gbk')
    parser.add_argument('gff_file', type=argparse.FileType("r"), help='GFF3 file')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), help='Fasta Input')
    args = parser.parse_args()

    for record in gff3_to_genbank(**vars(args)):
        SeqIO.write([record], sys.stdout, "genbank")
