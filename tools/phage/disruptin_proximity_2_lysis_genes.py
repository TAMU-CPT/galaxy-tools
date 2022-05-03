#!/usr/bin/env python
"""
This program is intended to identify protein coding sequences within a certain window (number of base pairs) of genes encoding recognized endolysin domains and genes encoding transmembrane domains. The goal is narrow a list of disruptin candidates by identifying the sequences close to the other lysis genes in the phage genome.
Inputs for this program include a .fasta file with protein sequences of lysis gene candidates from the phage genome, a .gff3 file with the tmhmm results from the genome, a .gff3 file with the results from interproscan of the genome, a .gff3 file of the genome, window size in number of base pairs, a tab separated list of endolysin domains, and optional names of output files.
The program outputs lists of lysis gene candidates that are close to protein codings sequences with endolysin domains or to sequences with transmembrane domains and lists of the proteins in proximity to the lysis gene candidates (one list for proteins with endolysin domains and one list for TMD-containing proteins).
"""

from Bio import SeqIO
import argparse
import sys
from cpt_gffParser import gffParse, gffWrite
from BCBio.GFF import GFFExaminer
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq
from intervaltree import IntervalTree, Interval
from gff3 import feature_lambda, feature_test_type

# Used for genome in fasta format
# outputs the start and end coordinates from a record in fasta format
# Should work for seqrecords from NCBI database
# def location_start_end(seqrec):
#    F = seqrec.description
#    description_subparts = F.split(' ')
#    for i in range(len(description_subparts)):
#        if description_subparts[i].startswith('[location'):
#            location = description_subparts[i][10:-1]
#    location_se = location.split('..')
#    location_start = location_se[0]
#    location_end = location_se[1]
#
#    return location_start, location_end

# adapted from intersect_and_adjacent.py
def treeFeatures(features):
    for feat in features:
        # used with genome in fasta format
        # start, end = location_start_end(feat)

        # Interval(begin, end, data)
        yield Interval(int(feat.location.start), int(feat.location.end), feat.id)


# Function to read enzyme domain names and ids from the enzyme list
# Enzyme list must be a tab separated txt file with the format with four columns: Predicted catalytic or binding domain, Abbreviation, Conserved domain, Phage example
# The first column is used here as the domain name, and the 3rd column is the domain id
def read_enzyme_list(enzyme_file=None):
    enzyme_file.seek(0)
    domains = []
    #domain_names = []

    for line in enzyme_file:
        if not line.startswith("*"):
            words = line.split("\t")
            if len(words) > 3:
                domains += [words[2]]
                #domain_names += [words[0]]

    return (domains[1:])


# adapted from intersect_and_adjacent.py
def intersect(rec_a, rec_b, window):
    if len(rec_a) > 0 and len(rec_b) > 0:

        # builds interval tree from Interval objects of form (start, end, id) for each feature
        tree_a = IntervalTree(list(treeFeatures(rec_a)))
        tree_b = IntervalTree(list(treeFeatures(rec_b)))

        # Used to map ids back to features later
        rec_a_map = {f.id: f for f in rec_a}
        rec_b_map = {f.id: f for f in rec_b}

        rec_a_hits_in_b = []
        rec_b_hits_in_a = []

        for feature in rec_a:
            # Used with genome in fasta format
            # start, end = location_start_end(feature)
            # Save each feature in rec_a that overlaps a feature in rec_b
            # hits = tree_b.find_range((int(feature.location.start), int(feature.location.end)))
            hits = tree_b[
                (int(feature.location.start) - window) : (
                    int(feature.location.end) + window
                )
            ]
            # feature id is saved in interval result.data, use map to get full feature
            for hit in hits:
                rec_a_hits_in_b.append(rec_b_map[hit.data])

        for feature in rec_b:
            # Used with genome in fasta format
            # start, end = location_start_end(feature)
            # Save each feature in rec_a that overlaps a feature in rec_b
            # hits = tree_a.find_range((int(feature.location.start), int(feature.location.end)))
            hits = tree_a[
                (int(feature.location.start) - window) : (
                    int(feature.location.end) + window
                )
            ]
            # feature id is saved in interval result.data, use map to get full feature
            for hit in hits:
                rec_b_hits_in_a.append(rec_a_map[hit.data])

        # Remove duplicate features using sets
        rec_a = set(rec_a_hits_in_b)
        rec_b = set(rec_b_hits_in_a)

    else:
        # If one input is empty, output two empty result files.
        rec_a = SeqRecord(Seq(""), "none")
        rec_b = SeqRecord(Seq(""), "none")
    return rec_a, rec_b


# Function to identify enzyme domains associated with endolysin function from the interproscan results file
def find_endolysins(rec_ipro, enzyme_domain_ids):

    # print(rec_ipro)

    if len(rec_ipro) > 0:
        endo_rec_names = []
        endo_rec_domain_ids = []
        rec_domain_name = []

        # Check each feature in the ipro file if the domain id/feature qualifier "Name" is included in the domain list.
        for seq in rec_ipro:
            for i in range(len(seq.features)):
                f = seq.features[i]
                if f.type == "protein_match":
                    # print(f)
                    unwanted = ["TRANS", "SIGNAL", "CYTO", "Coil", "Signal"]
                    # Ignores feature with unwanted key words in the feature name
                    if all(x not in f.qualifiers["Name"][0] for x in unwanted):
                        # If feature is included in the given enzyme domain list, the protein name, domain id, and domain name are stored
                        domain_description = [str(f.qualifiers['Name'][0])]
                        if 'signature_desc' in f.qualifiers:
                            domain_description += [str(f.qualifiers['signature_desc'][0])]

                        for i in enzyme_domain_ids:
                            for y in domain_description:
                                if i in y:
                                    endo_rec_domain_ids += [f.qualifiers["Name"][0]]

                                    # e_index = enzyme_domain_ids.index(i)
                                    # rec_domain_name += [enzyme_domain_names[e_index]]

                                    target = f.qualifiers["Target"][0]
                                    target = target.split(" ")
                                    protein_name = str(target[0]) + '**'
                                    endo_rec_names += [protein_name]

        return endo_rec_names, endo_rec_domain_ids


def adjacent_lgc(lgc, tmhmm, ipro, genome, enzyme, window):
    rec_lgc = list(SeqIO.parse(lgc, "fasta"))
    rec_tmhmm = gffParse(tmhmm)
    rec_ipro = gffParse(ipro)
    recTemp = gffParse(genome)[0]
    tempFeats = feature_lambda(
                recTemp.features,
                feature_test_type,
                {"types": ["CDS"]},
                subfeatures=True,
                )
    recTemp.features = tempFeats
    rec_genome_ini = [recTemp]

    # genome.seek(0)
    # examiner = GFFExaminer()
    # print(examiner.available_limits(genome))

    enzyme_domain_ids = read_enzyme_list(enzyme)

    if len(rec_lgc) > 0 and len(rec_tmhmm) > 0 and len(rec_genome_ini) > 0:

        # find names of the proteins containing endolysin associated domains
        endo_names, endo_domain_ids = find_endolysins(
            rec_ipro, list(enzyme_domain_ids)
        )

        # find names of proteins containing transmembrane domains
        tmhmm_protein_names = []
        for seq in rec_tmhmm:
            tmhmm_protein_names += [str(seq.id) + '**']

        lgc_names = []
        for seq in rec_lgc:
            lgc_names += [str(seq.id) + '**']

        adjacent_endo = {}
        adjacent_lgc_to_endo = {}
        adjacent_tm = {}
        adjacent_lgc_to_tm = {}

        # print(len(tmhmm_protein_names), len(endo_names))
        # print(rec_genome_ini)
        # print(len(rec_genome_ini))

        for i in range(len(rec_genome_ini)):
            rec_genome = rec_genome_ini[i]

            # find records for proteins containing endolysin domains and tmds from genome fasta file
            tm_seqrec = []
            endolysin_seqrec = []
            lgc_seqrec = []

            # print(tmhmm_protein_names)
            # print(endo_names)
            # print(lgc_names)
            # print(rec_genome)

            for feat in rec_genome.features:
                # print(feat)
                # searches for synonyms and
                if feat.type == "CDS":
                    feat_names = []
                    feat_names.append(str(feat.id) + '**')
                    if "locus_tag" in feat.qualifiers:
                        feat_names.append(str(feat.qualifiers["locus_tag"][0]) + '**')
                    if "protein_id" in feat.qualifiers:
                        feat_names.append(str(feat.qualifiers["protein_id"][0]) + '**')
                    if "Name" in feat.qualifiers:
                        if len(str(feat.qualifiers["Name"][0])) > 5:
                            feat_names.append(str(feat.qualifiers["Name"][0]) + '**')

                    # print(str(feat_names))
                    # print(str(feat.qualifiers))
                    for i in range(len(feat_names)):
                        if str(feat_names[i]) in str(lgc_names):
                            lgc_seqrec += [feat]
                    # check if gene annotated as holin using key words/synonyms
                    holin_annotations = ["holin"]
                    if "product" in feat.qualifiers:
                        if any(
                            x
                            for x in holin_annotations
                            if (x in str(feat.qualifiers))
                        ):
                            tm_seqrec += [feat]
                    # check if protein contains a TMD
                    for i in range(len(feat_names)):
                        if str(feat_names[i]) in str(tmhmm_protein_names):
                            # print(feat_names[i])
                            tm_seqrec += [feat]

                    # check if gene annotated as endolysin using key words/synonyms
                    endolysin_annotations = ["lysin", "lysozyme"]
                    if "product" in feat.qualifiers:
                        if any(
                            x
                            for x in endolysin_annotations
                            if (x in str(feat.qualifiers))
                        ):
                            endolysin_seqrec += [feat]
                    # check if protein contains an endolysin-associated domain
                    for i in range(len(feat_names)):
                        if str(feat_names[i]) in str(endo_names):
                            endolysin_seqrec += [feat]

            # print(endolysin_seqrec, tm_seqrec, lgc_seqrec)
            # find possible endolysins that are adjacent to (or within window length away from) the lysis gene, or disruptin, candidates
            # if len(endolysin_seqrec) > 0:
            adjacent_lgc_to_endo_i, adjacent_endo_i = intersect(
                endolysin_seqrec, lgc_seqrec, window
            )
            # find TMD-containing proteins that are adjacent to (or within window length away from) the lysis gene, or disruptin, candidates
            # if len(tm_seqrec) > 0:
            adjacent_lgc_to_tm_i, adjacent_tm_i = intersect(
                tm_seqrec, lgc_seqrec, window
            )

            # print(len(endolysin_seqrec), len(lgc_seqrec), len(tm_seqrec))
            adjacent_endo[rec_genome.id] = adjacent_endo_i
            adjacent_lgc_to_endo[rec_genome.id] = adjacent_lgc_to_endo_i
            adjacent_tm[rec_genome.id] = adjacent_tm_i
            adjacent_lgc_to_tm[rec_genome.id] = adjacent_lgc_to_tm_i
            # print(rec_genome.id)
    else:
      return 0, 0, 0, 0
    # print(adjacent_endo)
    return adjacent_endo, adjacent_lgc_to_endo, adjacent_tm, adjacent_lgc_to_tm


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Identify lysis gene candidates next to possible endolysin or holin genes",
        epilog="",
    )
    parser.add_argument(
        "lgc",
        type=argparse.FileType("r"),
        help="fasta file with protein coding sequences of lysis gene candidates",
    )
    parser.add_argument(
        "tmhmm",
        type=argparse.FileType("r"),
        help="gff3 file with tmhmm results from the coding sequences in the genome",
    )
    parser.add_argument(
        "ipro",
        type=argparse.FileType("r"),
        help="gff3 file with interpro results from protein sequences in the genome",
    )
    parser.add_argument(
        "genome",
        type=argparse.FileType("r"),
        help="fasta file with protein coding sequences for all genes in the genome",
    )
    parser.add_argument(
        "window",
        type=int,
        default=1000,
        help="Allows features this far away to still be considered 'adjacent'",
    )
    parser.add_argument(
        "enzyme",
        type=argparse.FileType("r"),
        help="tab delimited text file including list of conserved protein domains linked to endolysin function",
    )
    parser.add_argument("--oa", type=str, default="possible_endolysin.gff3")
    parser.add_argument(
        "--ob", type=str, default="lysis_gene_candidates_near_endolysin.gff3"
    )
    parser.add_argument("--oc", type=str, default="possible_holin.gff3")
    parser.add_argument(
        "--od", type=str, default="lysis_gene_candidates_near_holin.gff3"
    )
    args = parser.parse_args()

    endo, lgc_endo, tm, lgc_tm = adjacent_lgc(
        args.lgc, args.tmhmm, args.ipro, args.genome, args.enzyme, args.window
    )

    if endo == 0:
      with open(args.oa, "w") as handle:
        handle.write("##gff-version 3")
      with open(args.ob, "w") as handle:
        handle.write("##gff-version 3")
      with open(args.oc, "w") as handle:
        handle.write("##gff-version 3")
      with open(args.od, "w") as handle:
        handle.write("##gff-version 3")
      return

    args.genome.seek(0)
    rec = list(gffParse(args.genome))
    
    with open(args.oa, "w") as handle:
        for i in range(len(rec)):
            rec_i = rec[i]
            if endo.get(rec_i.id, "") is not "":
                rec_i.features = endo[rec_i.id]
                gffWrite([rec_i], handle)

    with open(args.ob, "w") as handle:
        for i in range(len(rec)):
            rec_i = rec[i]
            if lgc_endo.get(rec_i.id, "") is not "":
                rec_i.features = lgc_endo[rec_i.id]
                gffWrite([rec_i], handle)

    with open(args.oc, "w") as handle:
        for i in range(len(rec)):
            rec_i = rec[i]
            if tm.get(rec_i.id, "") is not "":
                rec_i.features = tm[rec_i.id]
                gffWrite([rec_i], handle)

    with open(args.od, "w") as handle:
        for i in range(len(rec)):
            rec_i = rec[i]
            if lgc_tm.get(rec_i.id, "") is not "":
                rec_i.features = lgc_tm[rec_i.id]
                gffWrite([rec_i], handle)
