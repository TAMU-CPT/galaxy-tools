#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
import logging
logging.basicConfig(level=logging.INFO)


def glimmer3_to_gff3(glimmer, genome):
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))

    current_record = None
    for line in glimmer:
        if line.startswith('>'):
            chromId = line.strip().replace('>', '')
            if chromId in seq_dict:
                if current_record is not None:
                    yield current_record
                current_record = seq_dict[chromId]
            else:
                raise Exception("Found results for sequence %s which was not in fasta file sequences (%s)" % (chromId, ', '.join(seq_dict.keys())))

        if not line.startswith('>'):
            (gene_id, gstart, gend, phase, score) = line.strip().split()
            gstart = int(gstart)
            gend = int(gend)

            if '+' in phase:
                strand = 1
                start = gstart
                end = gend
            else:
                strand = -1
                start = gend
                end = gstart

            # Correct for gff3
            start -= 1

            cds_feat = SeqFeature(
                FeatureLocation(start, end),
                type="CDS",
                strand=strand,
                qualifiers={
                    'source': 'Glimmer3',
                    'ID': '%s.cds_%s' % (current_record.id, gene_id),
                }
            )

            gene = SeqFeature(
                FeatureLocation(start, end),
                type="gene",
                strand=strand,
                qualifiers={
                    'source': 'glimmer',
                    'ID': '%s.%s' % (current_record.id, gene_id),
                }
            )
            gene.sub_features = [cds_feat]
            current_record.features.append(gene)
    yield current_record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Glimmer to GFF3')
    parser.add_argument('glimmer', type=argparse.FileType("r"), help='Glimmer3 Output')
    parser.add_argument('genome', type=argparse.FileType("r"), help='Fasta Genome')
    args = parser.parse_args()

    for result in glimmer3_to_gff3(**vars(args)):
        GFF.write([result], sys.stdout)
