#!/usr/bin/env python
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import logging
logging.basicConfig(level=logging.INFO)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), help='Fasta Input')
    parser.add_argument('sliceLen', type=int, help="Number of bases per slice")
    args = parser.parse_args()

    for sequence in SeqIO.parse(args.fasta_file, 'fasta'):
      j = 0
      for i in range(int(len(sequence.seq) / args.sliceLen)):
        temp_id = (sequence.id + "_" + str(i+1))
        temp_seq = (str(sequence.seq[i * args.sliceLen:((i + 1) * args.sliceLen)]))
        temp_rec = SeqRecord(Seq(temp_seq), id = temp_id, description = "")
        SeqIO.write(temp_rec, sys.stdout, 'fasta')
        print("")
        j += 1
      if (args.sliceLen == len(sequence.seq)):
        break
      temp_id = (sequence.id + "_" + str(j+1))
      temp_seq = (str(sequence.seq[j * args.sliceLen:]))
      temp_rec = SeqRecord(Seq(temp_seq), id = temp_id, description = "")
      SeqIO.write(temp_rec, sys.stdout, 'fasta')
      break # Only Single fasta for now
