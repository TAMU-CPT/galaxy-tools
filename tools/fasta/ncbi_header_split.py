#!/usr/bin/env python
import sys
import copy
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def split(fasta_file=None, id="merged"):
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        for header in rec.description.split('>'):
            nrec = copy.deepcopy(rec)
            nrec.id = header
            nrec.description = ''
            yield nrec


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', type=file)

    args = parser.parse_args()
    for seq in split(**vars(args)):
        SeqIO.write([seq], sys.stdout, 'fasta')
