#!/usr/bin/env python
import os
import sys
import argparse
from Bio import SeqIO
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split GBK file')
    parser.add_argument('genbank_file', type=file, help='Genbank file')
    args = parser.parse_args()

    outdir = os.path.join(os.getcwd(), 'gbk_out')
    os.makedirs(outdir)

    for record in SeqIO.parse(args.genbank_file, "genbank"):
        name = os.path.join(outdir, '%s.gbk' % record.id)
        if not os.path.exists(name):
            with open(name, 'w') as handle:
                log.info("Storing %s to %s", record.id, name)
                SeqIO.write([record], handle, 'genbank')
