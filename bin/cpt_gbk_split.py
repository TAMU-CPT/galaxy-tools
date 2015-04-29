#!/usr/bin/env python
import os
import argparse
from Bio import SeqIO
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split GBK file')
    parser.add_argument('genbank_file', type=file, help='Genbank file')
    args = parser.parse_args()

    for record in SeqIO.parse(args.genbank_file, "genbank"):
        name = '%s.gbk' % record.id
        if not os.path.exists(name):
            with open('%s.gbk' % record.id, 'w') as handle:
                log.info("Storing %s" % record.id)
                SeqIO.write([record], handle, 'genbank')
