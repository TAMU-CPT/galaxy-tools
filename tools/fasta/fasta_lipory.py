#!/usr/bin/env python
import re
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

LIPOBOX = re.compile('[ILMFTV][^REKD][GAS]C')


def get_id(feature=None, parent_prefix=None, idx=None):
    result = ""
    if parent_prefix is not None:
        result += parent_prefix + '|'
    if idx is not None:
        result += '%03d_' % idx
    if 'locus_tag' in feature.qualifiers:
        result += feature.qualifiers['locus_tag'][0]
    elif 'gene' in feature.qualifiers:
        result += feature.qualifiers['gene'][0]
    elif 'product' in feature.qualifiers:
        result += feature.qualifiers['product'][0]
    else:
        result += '%s_%s_%s' % (feature.location.start, feature.location.end,
                                feature.location.strand)
    return result


def find_lipoprotein(genbank_file=None, lipobox_mindist=10, lipobox_maxdist=30):
    for record in SeqIO.parse(genbank_file, "genbank"):
        if LIPOBOX.search(record.seq):

            newseq = ""
            last_idx = 0
            should_yield = False
            for hit in LIPOBOX.finditer(record.seq):
                if lipobox_mindist < hit.start() < lipobox_maxdist:
                    newseq += record.seq[last_idx:hit.start()].lower()
                    newseq += record.seq[hit.start():hit.end()].upper()
                    last_idx = hit.end()
                    should_yield = True

            newseq += record.seq[last_idx:].lower()
            # Update sequence

            if should_yield:
                yield [SeqRecord(Seq(newseq), id=get_id(feature),
                                    description="")]


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Filter out lipoproteins', epilog="")
    parser.add_argument('genbank_file', type=file, help='Genbank file')
    parser.add_argument('--lipobox_mindist', type=int, help='Minimum distance in codons to start of lipobox', default=10)
    parser.add_argument('--lipobox_maxdist', type=int, help='Maximum distance in codons to start of lipobox', default=33)

    args = parser.parse_args()

    args = vars(parser.parse_args())
    for seq in find_lipoprotein(**args):
        SeqIO.write(seq, sys.stdout, 'fasta')
