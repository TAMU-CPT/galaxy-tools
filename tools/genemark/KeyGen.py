#!/usr/bin/env python
import argparse

def keyGen(fasta_file):
  f = open(fasta_file,"w+")
  f.write("ACTAAGGTACGTACAGTCCGATGTTTTTTTTTTTTTTTTCCGATCGGCCATCGTGTAGGACATTCCAGCTTGAAGTACGT\n521146385")
  f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Workaround for making GeneMark key')
    parser.add_argument('fasta_file', type=str, help='Fasta file')
    args = parser.parse_args()

    keyGen(**vars(args))
