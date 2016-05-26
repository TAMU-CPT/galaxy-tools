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
    for rec in SeqIO.parse(genome, 'fasta')
        for line in glimmer:
            if not line.startswith('>'):
                (id, gstart, gend, phase, score) = line.strip().split()
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
                        'source': 'glimmer',
                        'ID': 'cds_%s' % id
                    }
                )

                gene = SeqFeature(
                    FeatureLocation(start, end),
                    type="gene",
                    strand=strand,
                    qualifiers={
                        'source': 'glimmer',
                        'ID': id,
                    }
                )
                gene.sub_features = [cds_feat]
                rec.features.append(gene)
        yield rec

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Glimmer to GFF3')
    parser.add_argument('glimmer', type=file, help='Glimmer3 Output')
    parser.add_argument('genome', type=file, help='Fasta Genome')
    args = parser.parse_args()

    for result in glimmer3_to_gff3(**vars(args)):
        GFF.write([result], sys.stdout)
