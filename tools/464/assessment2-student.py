#!/usr/bin/env python
import sys
import argparse
import random
import logging
from datetime import date
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
import os
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='verify against expectations')
    parser.add_argument('user_email', help='User email')
    args = parser.parse_args()

    SeqIO.write(generate_sequence(**vars(args)), sys.stdout, 'fasta')
