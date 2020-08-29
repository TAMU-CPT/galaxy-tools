#!/usr/bin/env python
import sys
import logging
import argparse
import copy
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
from Bio.SeqFeature import SeqFeature, FeatureLocation

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def add_exons(rec):
    for gene in feature_lambda(
        rec.features, feature_test_type, {"type": "gene"}, subfeatures=True
    ):
        clean_gene = copy.deepcopy(gene)
        exon_start = None
        exon_end = None
        cds_list = []
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
            if exon_end is None:
                exon_end = cds.location.end
            exon_start = min(exon_start, cds.location.start)
            exon_end = max(exon_end, cds.location.end)
            cds_list.append(cds)
        if cds_list:
            # we found a CDS to adopt
            new_exon = SeqFeature(
                location=FeatureLocation(exon_start, exon_end),
                type="exon",
                qualifiers={
                    "source": ["cpt.prepApollo"],
                    "ID": ["%s.exon" % clean_gene.qualifiers["ID"][0]],
                    "Parent": clean_gene.qualifiers["ID"],
                },
                sub_features=cds_list,
            )
            for cds in cds_list:
                #update parent qualifier for cdss
                cds.qualifiers["Parent"] = new_exon.qualifiers["ID"]
            # gene.sub_features.append(new_exon)
            # get all the other children of gene that AREN'T a CDS including the new exon
            clean_gene.sub_features = [copy.deepcopy(new_exon)]
            for sf in feature_lambda(
                gene.sub_features,
                feature_test_type,
                {"type": "CDS"},
                subfeatures=True,
                recurse=False,
                invert=True,
            ):
                child = copy.deepcopy(sf)
                child.qualifiers["Parent"] = clean_gene.qualifiers["ID"]
                clean_gene.sub_features.append(child)
            # add them to the gene feature
        # return the cleaned gene with new exon
        yield clean_gene


def gff_filter(gff3):
    for rec in GFF.parse(gff3):
        rec.features = sorted(list(add_exons(rec)), key=lambda x: x.location.start)
        rec.annotations = {}
        GFF.write([rec], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="add parent exon features to CDSs for Apollo"
    )
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    gff_filter(**vars(args))
