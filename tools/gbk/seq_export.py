#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import sys
import argparse
from Bio import SeqIO
import cpt

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export sequences from Genbank files')
    parser.add_argument('genbank_files', nargs='+', type=argparse.FileType("r"), help='Genbank file')
    parser.add_argument('--name_src', choices=('id', 'name', 'phage_name'), default='id')
    args = parser.parse_args()

    for gbk in args.genbank_files:
        for seq in SeqIO.parse(gbk, 'genbank'):
            if args.name_src == 'id':
                pass
            elif args.name_src == 'phage_name':
                (host, phage) = cpt.phage_name_parser(seq.description)
                seq.id = phage
            if args.name_src == 'name' or phage == None:
                seq.id = seq.name
            
                
            SeqIO.write(seq, sys.stdout, 'fasta')
