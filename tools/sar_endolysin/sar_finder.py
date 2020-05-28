#!/usr/bin/env python

import sys
import argparse
import os
import re
import pandas as pd
import numpy as np
from biopython_parsing import FASTA_parser
from Bio import SeqIO

class CheckSequence:
    """ SAR endolysin Verification class, which starts with complete FA file, and is shrunk by each function to reveal best candidates of SAR endolysin proteins """

    def __init__(self, protein_name, protein_data):
        self.name = protein_name
        self.seq = protein_data.seq
        self.size = len(self.seq)

    def check_sizes(self,min,max):
        """ check the minimum and maximum peptide lengths """
        if self.size < min:
            print("too small")
        elif self.size > max:
            print("too large")
        else:
            print(f"{self.name} : {self.seq}")
            return True

    def check_hydrophobicity(self):
        """ verifies the existence of a hydrophobic region within the sequence """
        hydrophobic_residues = "['FIWLVMYCATGS']"
        loc = 0
        tmd_size = 15 # current implementation is to have the TMD span 15 residues
        if self.size > 50:
            seq = self.seq[0:50]
        else:
            seq = self.seq
        for i in range(0,len(seq)-tmd_size,1):
            check_seq = str(seq[i:loc+tmd_size])
            print(f"original sequence: {seq}\ncheck sequence: {check_seq}")
            if re.search((hydrophobic_residues+"{"+str(tmd_size)+"}"),check_seq):
                print("hydrophobic region FOUND")
                return True
            else:
                print("did not find hydrophobic region")
                loc += 1
                continue


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SAR Finder")

    parser.add_argument("fa",type=argparse.FileType("r"),help="organism's multi fasta file")

    parser.add_argument("--min",type=int,default=20,help="minimum size of candidate peptide")

    parser.add_argument("--max",type=int,default=200,help="maximum size of candidate peptide")

    args = parser.parse_args()

    fa_dict = FASTA_parser(fa=args.fa).multifasta_dict()

    for protein_name, protein_data in fa_dict.items():
        sar = CheckSequence(protein_name, protein_data)
        #sar.check_sizes(min=args.min,max=args.max)
        sar.check_hydrophobicity()
