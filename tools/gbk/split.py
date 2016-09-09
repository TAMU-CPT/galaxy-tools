#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import os
import random
import argparse
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split GBK file')
    parser.add_argument('genbank_files', type=file, nargs='+', help='Genbank file')
    parser.add_argument('--allow_dupes', action='store_true', help='Allow files with duplicate IDs')
    args = parser.parse_args()

    outdir = os.path.join(os.getcwd(), 'gbk_out')
    os.makedirs(outdir)

    for genbank_file in args.genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            name = os.path.join(outdir, '%s.gbk' % record.id)

            # In case of dupes
            if os.path.exists(name) and not args.allow_dupes:
                continue

            while os.path.exists(name):
                random_id = ''.join([random.choice('1234567890') for x in range(12)])
                name = os.path.join(outdir, '%s.%s.gbk' % (record.id, random_id))

            with open(name, 'w') as handle:
                log.info("Storing %s to %s", record.id, name)
                SeqIO.write([record], handle, 'genbank')
