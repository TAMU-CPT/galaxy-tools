#!/usr/bin/env python
import sys
import logging
import argparse
import copy
from cpt_gffParser import gffParse, gffWrite, gffSeqFeature
from gff3 import feature_lambda, feature_test_type
from Bio.SeqFeature import FeatureLocation

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

ALLOWED_FEATURES = [
        "mRNA",
        "exon",
        "transposable_element",
        "tRNA",
        "transcript",
        "terminator",
        "Shine_Dalgarno_Sequence",
        "pseudogene",
        "stop_codon_read_through",
        "repeat_region",
        "CDS",
        "gene",
        "rRNA",
        "ncRNA",
        "snRNA",
        "snoRNA",
        "miRNA",
        ]

SPECIAL_REMOVED_FEATURES = ["gene_component_region", "sequence_difference"]



def add_exons(features):
    for gene in feature_lambda(
        features, feature_test_type, {"type": "gene"}, subfeatures=True
    ):
        clean_gene = copy.deepcopy(gene)
        exon_start = None
        exon_end = None
        exon_strand = None
        cds_list = []

        for mRNA in gene.sub_features:
            for x in mRNA.sub_features:
                x.qualifiers["Parent"] = [gene.ID]
                gene.sub_features.append(x)
                 
        for exon in feature_lambda(gene.sub_features, feature_test_type, {"type": "exon"}, subfeatures=False,recurse=False):
            #if the gene contains an exon, skip.
            continue
        
        # check for CDS child features of the gene, do not go a further step (this should skip any CDS children of exon child features)
        for cds in feature_lambda(
            gene.sub_features,
            feature_test_type,
            {"type": "CDS"},
            subfeatures=False,
            recurse=False,
        ):
            # check all CDS features for min/max boundaries
            if exon_start is None:
                exon_start = cds.location.start
                exon_strand = cds.location.strand
            if exon_end is None:
                exon_end = cds.location.end
            exon_start = min(exon_start, cds.location.start)
            exon_end = max(exon_end, cds.location.end)
            cds_list.append(cds)
        if cds_list:
            # we found a CDS to adopt
            new_exon = gffSeqFeature(
                location=FeatureLocation(exon_start, exon_end),
                type="exon",
                source = "cpt.prepApollo",
                qualifiers={
                    "ID": ["%s.exon" % clean_gene.qualifiers["ID"][0]],
                    "Parent": clean_gene.qualifiers["ID"],
                },
                sub_features=[],
                strand=exon_strand
            )
            for cds in cds_list:
                cds.qualifiers["Parent"] = new_exon.qualifiers["ID"]
            # gene.sub_features.append(new_exon)
            # get all the other children of gene that AREN'T a CDS including the new exon
            #clean_gene.sub_features = [copy.deepcopy(new_exon)]
            clean_gene.sub_features.append(gffSeqFeature(location=FeatureLocation(exon_start, exon_end, exon_strand), type="exon", source = "cpt.prepApollo", qualifiers={"ID": ["%s.exon" % clean_gene.qualifiers["ID"][0]], "Parent": clean_gene.qualifiers["ID"]}, sub_features=[], strand=exon_strand))
            """
            for sf in feature_lambda(
                gene.sub_features,
                feature_test_type,
                {"type": "CDS"},
                subfeatures=True,
                recurse=False,
                invert=True,
            ):
                child = copy.deepcopy(sf)
                child.qualifiers["Parent"] = new_exon.qualifiers["ID"]
                clean_gene.sub_features.append(child)
            """
            # add them to the new Exon feature
        # return the cleaned gene with new exon
        yield clean_gene

def process_features(features):
    # change RBS to 'Shine_Dalgarno_sequence'
    for rbs in feature_lambda(features, feature_test_type, {'type': "RBS"}):
        rbs.type = "Shine_Dalgarno_sequence"

    # Filter top level features
    for feature in feature_lambda(features, feature_test_type, {"types": ALLOWED_FEATURES}, subfeatures=True):
        cleaned_subfeatures = []
        for sf in feature.sub_features:
            if sf.type in SPECIAL_REMOVED_FEATURES:
                # 'gene_component_region' is uncaught by feature_test_type as it contains `gene`
                continue
            else:
                cleaned_subfeatures.append(sf)
        feature.sub_features = copy.deepcopy(cleaned_subfeatures)  
        yield feature

def gff_filter(gff3):
    for rec in gffParse(gff3):
        cleaned_features = sorted(list(process_features(rec.features)), key=lambda x: x.location.start)
        rec.features = sorted(list(add_exons(cleaned_features)), key=lambda x: x.location.start)
        rec.annotations = {}
        gffWrite([rec], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="add parent exon features to CDSs for Apollo"
    )
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    gff_filter(**vars(args))
