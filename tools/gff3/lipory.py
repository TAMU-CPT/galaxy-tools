#!/usr/bin/env python
import re
import sys
import argparse
import logging
from Bio import SeqIO
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type, get_id
from Bio.SeqFeature import SeqFeature, FeatureLocation
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def find_lipoprotein(gff3_file, fasta_genome, lipobox_mindist=10, lipobox_maxdist=60):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_genome, "fasta"))

    CASES = [
        re.compile('^.{%s,%s}[ILMFTV][^REKD][GAS]C' % (lipobox_mindist, lipobox_maxdist)),
        re.compile('^.{%s,%s}AWAC' % (lipobox_mindist, lipobox_maxdist)),
    ]

    for record in GFF.parse(gff3_file, base_dict=seq_dict):
        good_features = []

        genes = list(feature_lambda(record.features, feature_test_type, {'type': 'gene'}, subfeatures=True))
        for gene in genes:
            cdss = list(feature_lambda(gene.sub_features, feature_test_type, {'type': 'CDS'}, subfeatures=False))
            if len(cdss) == 0:
                continue

            # Someday this will bite me in the arse.
            cds = cdss[0]

            tmpseq = str(cds.extract(record.seq).translate(table=11, cds=True)).replace("*", "")
            for case in CASES:
                m = case.search(tmpseq)
                if m:
                    if cds.location.strand > 0:
                        start = cds.location.start + m.end() - 4
                        end = cds.location.start + m.end()
                    else:
                        start = cds.location.end - m.end() + 4
                        end = cds.location.end - m.end()

                    tmp = SeqFeature(
                        FeatureLocation(
                            min(start, end),
                            max(start, end),
                            strand=cds.location.strand
                        ),
                        type='Lipobox',
                        qualifiers={
                            'source': 'CPT_LipoRy',
                            'ID': '%s.lipobox' % get_id(gene)
                        }
                    )

                    gene.sub_features.append(tmp)
                    good_features.append(gene)

        record.features = good_features
        yield [record]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter out lipoproteins', epilog="")
    parser.add_argument('gff3_file', type=file, help='Naive ORF Calls')
    parser.add_argument('fasta_genome', type=file, help='Fasta genome sequence')

    parser.add_argument('--lipobox_mindist', type=int,
                        help='Minimum distance in codons to start of lipobox', default=10)
    parser.add_argument('--lipobox_maxdist', type=int,
                        help='Maximum distance in codons to start of lipobox', default=33)

    args = parser.parse_args()

    args = vars(parser.parse_args())
    for record in find_lipoprotein(**args):
        record[0].annotations = {}
        GFF.write(record, sys.stdout)
