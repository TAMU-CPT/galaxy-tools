#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
import logging
logging.basicConfig(level=logging.INFO)


def mga_to_gff3(mga_output, genome):
    rec = SeqIO.read(genome, 'fasta')

    for line in mga_output:
        if not line.startswith('#'):
            (gene_id, start, end, strand, phase, complete, score, model,
                rbs_start, rbs_end, rbs_score) = line.strip().split('\t')
            start = int(start)
            end = int(end)
            strand = +1 if strand == '+' else -1

            # Correct for gff3
            start -= 1

            rbs_feat = None
            if rbs_start != '-':
                rbs_start = int(rbs_start)
                rbs_end = int(rbs_end)
                rbs_feat = SeqFeature(
                    FeatureLocation(rbs_start, rbs_end),
                    type="Shine_Dalgarno_sequence",
                    strand=strand,
                    qualifiers={
                        'ID': 'rbs_%s' % gene_id
                    }
                )

            cds_feat = SeqFeature(
                FeatureLocation(start, end),
                type="CDS",
                strand=strand,
                qualifiers={
                    'ID': 'cds_%s' % gene_id
                }
            )

            if rbs_feat is not None:
                if strand > 0:
                    gene_start = rbs_start
                    gene_end = end
                else:
                    gene_start = start
                    gene_end = rbs_end
            else:
                gene_start = start
                gene_end = end

            gene = SeqFeature(
                FeatureLocation(gene_start, gene_end),
                type="gene",
                strand=strand,
                qualifiers={
                    'ID': gene_id
                }
            )

            gene.sub_features = [cds_feat]
            if rbs_feat is not None:
                gene.sub_features.append(rbs_feat)

            rec.features.append(gene)
    return rec


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert MGA to GFF3', epilog="")
    parser.add_argument('mga_output', type=file, help='MetaGeneAnnotator Output')
    parser.add_argument('genome', type=file, help='Fasta Genome')
    args = parser.parse_args()

    result = mga_to_gff3(**vars(args))

    result.annotations = {}
    GFF.write([result], sys.stdout)
