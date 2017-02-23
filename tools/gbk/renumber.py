#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
import sys  # noqa
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def renumber_genes(gbk_files, tag_to_update="locus_tag",
                   string_prefix="display_id", leading_zeros=3,
                   change_table=None):

    for gbk_file in gbk_files:
        for record in SeqIO.parse(gbk_file, "genbank"):
            if string_prefix == 'display_id':
                format_string = record.id + '_%0' + str(leading_zeros) + 'd'
            else:
                format_string = string_prefix + '%0' + str(leading_zeros) + 'd'

            f_cds = [f for f in record.features if f.type == 'CDS']
            f_rbs = [f for f in record.features if f.type == 'RBS']
            f_gene = [f for f in record.features if f.type == 'gene']
            f_oth = [f for f in record.features if f.type not in ['CDS', 'RBS',
                                                                  'gene']]
            f_care_about = []

            # Make sure we've hit every RBS and gene
            for cds in f_cds:
                # If there's an associated gene feature, it will share a stop codon
                if cds.location.strand > 0:
                    associated_genes = [f for f in f_gene if f.location.end ==
                                        cds.location.end]
                else:
                    associated_genes = [f for f in f_gene if f.location.start ==
                                        cds.location.start]

                # If there's an RBS it'll be upstream a bit.
                if cds.location.strand > 0:
                    associated_rbss = [f for f in f_rbs if f.location.end <
                                       cds.location.start and f.location.end >
                                       cds.location.start - 24]
                else:
                    associated_rbss = [f for f in f_rbs if f.location.start >
                                       cds.location.end and f.location.start <
                                       cds.location.end + 24]
                tmp_result = [cds]
                if len(associated_genes) > 0:
                    tmp_result.append(associated_genes[0])

                if len(associated_rbss) == 1:
                    tmp_result.append(associated_rbss[0])
                else:
                    log.warning("%s RBSs found for %s", len(associated_rbss), cds.location)
                # We choose to append to f_other as that has all features not
                # already accessed. It may mean that some gene/RBS features are
                # missed if they aren't detected here, which we'll need to handle.
                f_care_about.append(tmp_result)

            # Why must these people start at 1
            # Because we have to modify the array we work on a separate one
            clean_features = f_oth
            delta = []
            for index, feature_list in enumerate(sorted(f_care_about, key=lambda x: x[0].location.start)):
                for f in feature_list:
                    original_tag_value = delta_old(f, tag_to_update)
                    new_tag_value = format_string % index
                    f.qualifiers[tag_to_update] = [new_tag_value]
                    clean_features.append(f)
                    delta.append('\t'.join((record.id, original_tag_value, new_tag_value)))

            # Update all features
            record.features = sorted(clean_features, key=lambda x: x.location.start)

            change_table.write('\n'.join(delta) + '\n')

            # Output
            yield record


def delta_old(feature, tag_to_update):
    # First part of delta entry, old name
    if tag_to_update in feature.qualifiers:
        return feature.qualifiers[tag_to_update][0]
    else:
        return '%s %s %s' % (feature.location.start, feature.location.end,
                             feature.location.strand,)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Renumber genbank files')
    parser.add_argument('gbk_files', type=argparse.FileType("r"), nargs='+', help='Genbank files')
    parser.add_argument('--tag_to_update', type=str, help='Tag to update', default='locus_tag')
    parser.add_argument('--string_prefix', type=str, help='Prefix string', default='display_id')
    parser.add_argument('--leading_zeros', type=int, help='# of leading zeroes', default=3)

    parser.add_argument('--change_table', type=argparse.FileType('w'),
                        help='Location to store change table in', default='renumber.tsv')

    args = parser.parse_args()
    for record in renumber_genes(**vars(args)):
        SeqIO.write(record, sys.stdout, 'genbank')
