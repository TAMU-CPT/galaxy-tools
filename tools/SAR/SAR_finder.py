import sys
import argparse
import os
import re
from biopython_parsing import FASTA_parser
from file_operations import fasta_from_SAR_dict, gff3_from_SAR_dict, tab_from_SAR_dict
from SAR_functions import CheckSequence

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SAR Finder")

    parser.add_argument("fa",type=argparse.FileType("r"),help="organism's multi fasta file")

    parser.add_argument("--min",type=int,default=20,help="minimum size of candidate peptide")

    parser.add_argument("--max",type=int,default=200,help="maximum size of candidate peptide")

    parser.add_argument("--sar_min",type=int,default=15,help="minimum size of candidate peptide TMD domain")

    parser.add_argument("--sar_max",type=int,default=20,help="maximum size of candidate peptide TMD domain")
    
    parser.add_argument("--out_fa",type=argparse.FileType("w"),help="multifasta output of candidate SAR proteins",default="candidate_SAR.fa")

    parser.add_argument("--out_stat",type=argparse.FileType("w"),help="summary statistic file for candidate SAR proteins, tab separated",default="candidate_SAR_stats.tsv")

    parser.add_argument("--out_gff3",type=argparse.FileType("w"),help="multigff3 file for candidate SAR proteins",default="candidate_SAR.gff3")

    args = parser.parse_args()

    fa_dict = FASTA_parser(fa=args.fa).multifasta_dict()

    sars = {}
    for protein_name, protein_data in fa_dict.items():
        sar = CheckSequence(protein_name, protein_data)
        #sar.check_sizes(min=args.min,max=args.max)
        hydros = sar.shrink_results(sar_min=args.sar_min, sar_max=args.sar_max)
        sars.update(hydros)
    
    gff3_from_SAR_dict(sars, args.out_gff3)
    tab_from_SAR_dict(sars,args.out_stat,"SGA",sar_min=args.sar_min, sar_max=args.sar_max)
    fasta_from_SAR_dict(sars,args.out_fa)
    #stat_file_from_SAR_dict(sars,args.out_stat,sar_min=args.sar_min,sar_max=args.sar_max) # fix this whenever ready.