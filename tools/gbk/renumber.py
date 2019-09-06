#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
import sys  # noqa
from Bio import SeqIO

import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

# gene and RBS features are also included in the tagged features list, but are dealt with specifically elsewhere.
# This is used to filter out just valid "in gene" features
TAGGED_FEATURES = ["CDS", "tRNA", "intron", "mat_peptide"]


def renumber_genes(
    gbk_files,
    tag_to_update="locus_tag",
    string_prefix="display_id",
    leading_zeros=3,
    change_table=None,
):

    for gbk_file in gbk_files:
        for record in SeqIO.parse(gbk_file, "genbank"):
            if string_prefix == "display_id":
                format_string = record.id + "_%0" + str(leading_zeros) + "d"
            else:
                format_string = string_prefix + "%0" + str(leading_zeros) + "d"

            # f_cds = [f for f in record.features if f.type == 'CDS']
            # f_rbs = [f for f in record.features if f.type == 'RBS']
            # f_gene = [f for f in record.features if f.type == 'gene']
            # f_intron = [f for f in record.features if f.type == 'intron']
            # f_trna = [f for f in record.features if f.type == 'tRNA']
            # f_pep = [f for f in record.features if f.type == 'mat_peptide']
            # f_oth = [f for f in record.features if f.type not in ['CDS', 'RBS',
            #                                                      'gene', 'intron',
            #                                                      'tRNA', 'mat_peptide']]
            # Apparently we're numbering tRNAs now, thanks for telling me.
            # f_oth2 = []
            # for q in sorted(f_oth, key=lambda x: x.location.start):
            #    if q.type == 'tRNA':
            #        q.qualifiers['locus_tag'] = format_string_t % tRNA_count
            #        tRNA_count += 1
            #        f_oth2.append(q)
            #    else:
            #        f_oth2.append(q)
            # f_oth = f_oth2

            # f_care_about = []

            # Make sure we've hit every RBS and gene
            # for cds in f_cds:
            # If there's an associated gene feature, it will share a stop codon
            #    if cds.location.strand > 0:
            #        associated_genes = [f for f in f_gene if f.location.end ==
            #                            cds.location.end]
            #    else:
            #        associated_genes = [f for f in f_gene if f.location.start ==
            #                            cds.location.start]

            #    # If there's an RBS it'll be upstream a bit.
            #    if cds.location.strand > 0:
            #        associated_rbss = [f for f in f_rbs if f.location.end <
            #                           cds.location.start and f.location.end >
            #                           cds.location.start - 24]
            #    else:
            #        associated_rbss = [f for f in f_rbs if f.location.start >
            #                           cds.location.end and f.location.start <
            #                           cds.location.end + 24]
            #    tmp_result = [cds]
            #    if len(associated_genes) > 0:
            #        tmp_result.append(associated_genes[0])

            #   if len(associated_rbss) == 1:
            #       tmp_result.append(associated_rbss[0])
            #   else:
            #       log.warning("%s RBSs found for %s", len(associated_rbss), cds.location)
            # We choose to append to f_other as that has all features not
            # already accessed. It may mean that some gene/RBS features are
            # missed if they aren't detected here, which we'll need to handle.
            #    f_care_about.append(tmp_result)

            #####-----------------------------------------------------------------------------------------------------------#####
            # Build list of genes, then iterate over non-gene features and sort into containing genes.
            # tags are assigned based on genes, so start the lists with the gene features
            f_gene = sorted(
                [f for f in record.features if f.type == "gene"],
                key=lambda x: x.location.start,
            )
            f_rbs = sorted(
                [f for f in record.features if f.type == "RBS"],
                key=lambda x: x.location.start,
            )
            f_tag = list()
            f_sorted = sorted(
                [f for f in record.features if f.type in TAGGED_FEATURES],
                key=lambda x: x.location.start,
            )
            # genes not in the TAGGED_FEATURES list are exluded from the processing and assumed to already be clean
            clean_features = sorted(
                [
                    f
                    for f in record.features
                    if f.type not in TAGGED_FEATURES and f.type not in ["gene", "RBS"]
                ],
                key=lambda x: x.location.start,
            )

            f_processed = []
            for gene in f_gene:
                tag = [gene]
                # find the gene's RBS feature
                for rbs in [f for f in f_rbs if f not in f_processed]:
                    if is_within(rbs, gene):
                        tag.append(rbs)
                        f_processed.append(rbs)
                        break
                # find all other non-RBS features
                for feature in [f for f in f_sorted if f not in f_processed]:
                    # If the feature is within the gene boundaries (genes are the first entry in tag list),
                    # add it to the same locus tag group, does not process RBS
                    if is_within(feature, gene):
                        # catches genes and CDS feature that are intron-contained.
                        if feature.type == "CDS":
                            if (
                                feature.location.start == gene.location.start
                                or feature.location.end == gene.location.end
                            ):
                                tag.append(feature)
                                f_processed.append(feature)
                        else:
                            tag.append(feature)
                            f_processed.append(feature)
                    elif feature.location.start > gene.location.end:
                        # because the features are sorted by coordinates,
                        # no features further down  on the list will be in this gene
                        break
                f_tag.append(tag)

            # Process for frameshifts and mat_peptides (inteins)

            # check for overlapped genes
            # at this point, relevant features are put into tag buckets along with the containing gene
            # matin the form of [gene, feature1, feature2, ...]
            tag_index = 1
            delta = []
            for tag in f_tag:  # each tag list is one 'bucket'
                new_tag_value = format_string % tag_index
                for feature in tag:
                    original_tag_value = delta_old(feature, tag_to_update)
                    feature.qualifiers[tag_to_update] = [new_tag_value]
                    # Once the tag is renumbered, it's added to the clean list for later output
                    clean_features.append(feature)
                    delta.append(
                        "\t".join((record.id, original_tag_value, new_tag_value))
                    )
                tag_index += 1

            # Why must these people start at 1
            # Because we have to modify the array we work on a separate one
            # clean_features = f_oth
            # delta = []
            # for index, feature_list in enumerate(sorted(f_care_about, key=lambda x: x[0].location.start)):
            #    for f in feature_list:
            #        original_tag_value = delta_old(f, tag_to_update)
            #        # Add 1 to index for 1-indexed counting for happy scientists
            #        new_tag_value = format_string % (index+1)
            #        f.qualifiers[tag_to_update] = [new_tag_value]
            #        clean_features.append(f)
            #        delta.append('\t'.join((record.id, original_tag_value, new_tag_value)))

            # Update all features
            record.features = sorted(clean_features, key=lambda x: x.location.start)

            change_table.write("\n".join(delta) + "\n")

            # Output
            yield record


def delta_old(feature, tag_to_update):
    # First part of delta entry, old name
    if tag_to_update in feature.qualifiers:
        return feature.qualifiers[tag_to_update][0]
    else:
        return "%s %s %s" % (
            feature.location.start,
            feature.location.end,
            feature.location.strand,
        )


def is_within(query, feature):
    # checks if the query item is within the bounds of the given feature
    if (
        feature.location.start <= query.location.start
        and feature.location.end >= query.location.end
    ):
        return True
    else:
        return False


# def fix_frameshift(a, b):
#    #checks if gene a and gene b are a frameshifted gene (either shares a start or an end and an RBS)
#    if a[0].location.start == b[0].location.start or a[0].location.end == b[0].location.end:
#        # It is likely a frameshift. Treat is as such. Find shared RBS, determine which CDS is which
#        big_gene = a if (a[0].location.end - a[0].location.start) > (b[0].location.end - b[0].location.start) else b
#        small_gene = a if big_gene==b else b
#        rbs = [f for f in a if f.type == 'RBS']
#        # In the way that the tag lists are generated, the larger gene should contain both CDS features.
#        # Retrieve and dermine big/small CDS
#        cdss = [f for f in big_gene if f.type == 'CDS']
#        big_cds = cdss[0] if (cdss[0].location.end - cdss[0].location.start) > (cdss[1].location.end - cdss[1].location.start) else cdss[1]
#        small_cds = cdss[0] if big_cds==cdss[1] else cdss[1]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Renumber genbank files")
    parser.add_argument(
        "gbk_files", type=argparse.FileType("r"), nargs="+", help="Genbank files"
    )
    parser.add_argument(
        "--tag_to_update", type=str, help="Tag to update", default="locus_tag"
    )
    parser.add_argument(
        "--string_prefix", type=str, help="Prefix string", default="display_id"
    )
    parser.add_argument(
        "--leading_zeros", type=int, help="# of leading zeroes", default=3
    )

    parser.add_argument(
        "--change_table",
        type=argparse.FileType("w"),
        help="Location to store change table in",
        default="renumber.tsv",
    )

    args = parser.parse_args()
    for record in renumber_genes(**vars(args)):
        SeqIO.write(record, sys.stdout, "genbank")
