#!/usr/bin/env python
import sys
import argparse
from CPT_GFFParser import gffParse, gffWrite, gffSeqFeature
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
import logging

logging.basicConfig(level=logging.INFO)


def glimmer3_to_gff3(glimmer, genome):
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))

    current_record = None
    for line in glimmer:
        if line.startswith(">"):
            chromId = line.strip().replace(">", "")
            if chromId in seq_dict:
                if current_record is not None:
                    yield current_record
                current_record = seq_dict[chromId]
            else:
                raise Exception(
                    "Found results for sequence %s which was not in fasta file sequences (%s)"
                    % (chromId, ", ".join(seq_dict.keys()))
                )

        if not line.startswith(">"):
            (gene_id, gstart, gend, phase, score) = line.strip().split()
            gstart = int(gstart)
            gend = int(gend)

            if "+" in phase:
                strand = 1
                start = gstart
                end = gend
            else:
                strand = -1
                start = gend
                end = gstart

            # Correct for gff3
            start -= 1

            if start > end:
                #gene found on boundary (ex [4000, 200]) from glimmer assuming circular genome
                #-------------start<=======|sequence end|========>end------
                if strand > 0:
                    end = len(current_record)
                else:
                    start = 0
                gene_id+="_truncated"

            cds_feat = gffSeqFeature(
                FeatureLocation(start, end),
                type="CDS",
                strand=strand,
                id="%s.%s" % (current_record.id, gene_id),
                qualifiers={
                    "source": "Glimmer3",
                    "ID": "%s.cds_%s" % (current_record.id, gene_id),
                },
                source="Glimmer3"
            )

            gene = gffSeqFeature(
                FeatureLocation(start, end),
                type="gene",
                strand=strand,
                id="%s.%s" % (current_record.id, gene_id),
                qualifiers={
                    "source": "Glimmer3",
                    "ID": "%s.%s" % (current_record.id, gene_id),
                },
                source="Glimmer3"
            )
            gene.sub_features = [cds_feat]
            current_record.features.append(gene)
    yield current_record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Glimmer to GFF3")
    parser.add_argument("glimmer", type=argparse.FileType("r"), help="Glimmer3 Output")
    parser.add_argument("genome", type=argparse.FileType("r"), help="Fasta Genome")
    args = parser.parse_args()

    for result in glimmer3_to_gff3(**vars(args)):
        gffWrite([result], sys.stdout)
