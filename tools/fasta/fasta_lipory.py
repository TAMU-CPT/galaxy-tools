#!/usr/bin/env python
import re
import sys
import argparse
from Bio import SeqIO
from BCBio import GFF

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

LIPOBOX = re.compile('[ILMFTV][^REKD][GAS]C')


def find_lipoprotein(gff3_file, fasta_genome, lipobox_mindist=10, lipobox_maxdist=30):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_genome, "fasta"))

    for record in GFF.parse(gff3_file, base_dict=seq_dict):
        good_features = []
        for feature in record.features:
            if feature.type == 'remark':
                continue

            tmpseq = str(feature.extract(record.seq).translate(table=11)).replace("*", "")
            if LIPOBOX.search(tmpseq):
                good_features.append(feature)
        record.features = good_features
        yield [record]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter out lipoproteins', epilog="")
    parser.add_argument('gff3_file', type=file, help='Naive ORF Calls')
    parser.add_argument('fasta_genome', type=file, help='Fasta genome sequence')

    parser.add_argument('--lipobox_mindist', type=int, help='Minimum distance in codons to start of lipobox', default=10)
    parser.add_argument('--lipobox_maxdist', type=int, help='Maximum distance in codons to start of lipobox', default=33)

    args = parser.parse_args()

    args = vars(parser.parse_args())
    for record in find_lipoprotein(**args):
        GFF.write(record, sys.stdout)
