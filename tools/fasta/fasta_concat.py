#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def merge(fasta_file=None, id="merged"):
    sequence = ''
    ids = []
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for i in range(len(records)):
        ids.append(records[i].id)
        sequence += records[i].seq

    output = []
    output.append(SeqRecord(seq=sequence, id=id, description='Created from [%s...]' % (','.join(ids[0:10]))))
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), help='Genbank file')
    parser.add_argument('--id', help='New fasta identifier for merged sequences', default='merged')

    args = parser.parse_args()
    SeqIO.write(merge(**vars(args)), sys.stdout, "fasta")
