#!/usr/bin/env python
import logging
import sys
import argparse
from gff3 import feature_lambda
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import reverse_complement, translate
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


def feature_test(feature, **kwargs):
    return feature.type == kwargs['type'] and str(feature.extract(kwargs['record'].seq))[-3:] in kwargs['stops']


def suppress(genome, annotations, suppress=None):
    if suppress is None:
        raise Exception("Must provide a list of stop codons to suppress")

    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))

    suppressed_features = []

    for record in GFF.parse(annotations, base_dict=seq_dict):
        for feature in feature_lambda(
            record.features,
            feature_test,
            {'type': 'CDS', 'record': record, 'stops': suppress},
            subfeatures=True
        ):
            log.info("Found matching feature %s", feature.id)
            new_end = None

            codon_idx = 0
            while new_end is None:
                if feature.strand > 0:
                    cs = feature.location.end + (3 * codon_idx)
                    codon = str(record.seq[cs:cs + 3])
                else:
                    cs = feature.location.start - (3 * (1 + codon_idx))
                    codon = reverse_complement(record.seq[cs:cs + 3])

                if codon not in suppress and translate(codon, 11) == '*':
                    new_end = codon_idx
                    break

                codon_idx += 1
                if codon_idx > 40:
                    log.warn("Could not find a new stop codon")
                    break

            if new_end is not None:
                if feature.strand > 0:
                    feature.location._end += (codon_idx * 3)
                else:
                    feature.location._start -= (codon_idx * 3)
            suppressed_features.append(feature)

        record.features = suppressed_features
        record.annotations = {}
        GFF.write([record], sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='generates a genome with specified stop codons suppressed')
    parser.add_argument('genome', type=file, help='Genome Sequence')
    parser.add_argument('annotations', type=file, help='Parent GFF3 annotations')
    parser.add_argument('--suppress', type=str, nargs='+',
                        help=(
                            'Suppress this stop codon. All features with this '
                            'codon will be extended to the next available stop '
                            'codon'
                        ))
    args = parser.parse_args()

    suppress(**vars(args))
