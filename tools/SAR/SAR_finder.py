import sys
import argparse
import os
import re
import pandas as pd
import numpy as np
from biopython_parsing import FASTA_parser
from SAR_functions import CheckSequence
from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SAR Finder")

    parser.add_argument("fa",type=argparse.FileType("r"),help="organism's multi fasta file")

    parser.add_argument("--min",type=int,default=20,help="minimum size of candidate peptide")

    parser.add_argument("--max",type=int,default=200,help="maximum size of candidate peptide")

    args = parser.parse_args()

    fa_dict = FASTA_parser(fa=args.fa).multifasta_dict()

    hits = {}
    for protein_name, protein_data in fa_dict.items():
        sar = CheckSequence(protein_name, protein_data)
        #sar.check_sizes(min=args.min,max=args.max)
        hydros = sar.check_hydrophobicity()
        hits.update(hydros)
    
    print(hits)