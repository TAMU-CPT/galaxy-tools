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
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))

    current_record = None
    for line in mga_output:
        if line.startswith('#'):
            if line.startswith('# gc = ') or line.startswith('# self:'):
                continue
            chromId = line.strip().replace('# ', '')

            if ' ' in chromId:
                chromId = chromId[0:chromId.index(' ')]

            if chromId in seq_dict:
                if current_record is not None:
                    yield current_record
                current_record = seq_dict[chromId]
            else:
                raise Exception("Found results for sequence %s which was not in fasta file sequences (%s)" % (chromId, ', '.join(seq_dict.keys())))

        else:
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
                        'ID': '%s.rbs_%s' % (current_record.id, gene_id),
                        'Source': 'MGA'
                    }
                )

            cds_feat = SeqFeature(
                FeatureLocation(start, end),
                type="CDS",
                strand=strand,
                qualifiers={
                    'Source': 'MGA',
                    'ID': '%s.cds_%s' % (current_record.id, gene_id),
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
                    'Source': 'MGA',
                    'ID': '%s.%s' % (current_record.id, gene_id),
                }
            )

            gene.sub_features = [cds_feat]
            if rbs_feat is not None:
                gene.sub_features.append(rbs_feat)
            current_record.features.append(gene)
    yield current_record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert MGA to GFF3', epilog="")
    parser.add_argument('mga_output', type=argparse.FileType("r"), help='MetaGeneAnnotator Output')
    parser.add_argument('genome', type=argparse.FileType("r"), help='Fasta Genome')
    args = parser.parse_args()

    for result in mga_to_gff3(**vars(args)):
        GFF.write([result], sys.stdout)
