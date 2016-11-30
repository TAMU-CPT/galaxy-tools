#!/usr/bin/env python
import sys
import logging
import argparse
from Bio import SeqIO
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
from shinefind import NaiveSDCaller
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def require_shinefind(gff3, fasta):
    sd_finder = NaiveSDCaller()
    # Load up sequence(s) for GFF3 data
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    # Parse GFF3 records
    for record in GFF.parse(gff3, base_dict=seq_dict):
        # Reopen
        genes = list(feature_lambda(record.features, feature_test_type, {'type': 'gene'}, subfeatures=True))
        good_genes = []
        for gene in genes:
            cdss = list(feature_lambda(gene.sub_features, feature_test_type, {'type': 'CDS'}, subfeatures=False))
            if len(cdss) == 0:
                continue

            # Someday this will bite me in the arse.
            cds = cdss[0]

            sds, start, end, seq = sd_finder.testFeatureUpstream(cds, record, sd_min=5, sd_max=17)
            if len(sds) > 1:
                # TODO
                # Double plus yuck
                sd_features = sd_finder.to_features(sds, gene.location.strand, start, end, feature_id=gene.id)
                gene.sub_features.append(
                    sd_features[0]
                )

                good_genes.append(gene)

        # Yuck!
        record.features = good_genes
        yield record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta', type=file, help='Fasta Genome')
    parser.add_argument('gff3', type=file, help='GFF3 annotations')
    args = parser.parse_args()

    for rec in require_shinefind(**vars(args)):
        rec.annotations = {}
        GFF.write([rec], sys.stdout)
