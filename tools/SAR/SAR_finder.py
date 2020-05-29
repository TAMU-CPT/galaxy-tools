import sys
import argparse
import os
import re
from biopython_parsing import FASTA_parser
from file_operations import fasta_from_SAR_dict, stat_file_from_SAR_dict
from SAR_functions import CheckSequence

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SAR Finder")

    parser.add_argument("fa",type=argparse.FileType("r"),help="organism's multi fasta file")

    parser.add_argument("--min",type=int,default=20,help="minimum size of candidate peptide")

    parser.add_argument("--max",type=int,default=200,help="maximum size of candidate peptide")

    parser.add_argument("--tmd_min",type=int,default=15,help="minimum size of candidate peptide TMD domain")

    parser.add_argument("--tmd_max",type=int,default=20,help="maximum size of candidate peptide TMD domain")
    
    parser.add_argument("--out_fa",type=argparse.FileType("w"),help="multifasta output of candidate SAR proteins",default="candidate_SAR.fa")

    parser.add_argument("--out_stat",type=argparse.FileType("w"),help="summary statistic file for candidate SAR proteins",default="candidate_SAR_stats.txt")

    args = parser.parse_args()

    fa_dict = FASTA_parser(fa=args.fa).multifasta_dict()

    hits = {}
    for protein_name, protein_data in fa_dict.items():
        sar = CheckSequence(protein_name, protein_data)
        #sar.check_sizes(min=args.min,max=args.max)
        hydros = sar.check_hydrophobicity_and_charge(tmd_min=args.tmd_min, tmd_max=args.tmd_max)
        hits.update(hydros)
    
    fasta_from_SAR_dict(hits,args.out_fa)
    stat_file_from_SAR_dict(hits,args.out_stat,tmd_min=args.tmd_min,tmd_max=args.tmd_max)