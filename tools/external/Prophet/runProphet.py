#!/usr/bin/env python

import os
import argparse
import subprocess


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Python shell to suppress pipe errors in Prophet")

    parser.add_argument("--evalue_cutoff", type=str)
    parser.add_argument("--window_size", type=str)
    parser.add_argument("fasIn", type=str) 
    parser.add_argument("gffIn", type=str)

    args = parser.parse_args()
    cmd = "/galaxy/tools/cpt2/galaxy-tools/tools/external/Prophet/ProphET_standalone.pl --evalue_cutoff " + args.evalue_cutoff + " --window " + args.window_size + " --fasta " + args.fasIn + " --gff_in " + args.gffIn + " --outdir ./tempOut"
    try :
        subprocess.run(cmd, shell=True) # Just to keep any sub-module errors Prophet throws to keep Galaxy from thinking the job failed
    except :                            # The except is irrelevant, just needs to complete
        res = 0
    
    
    

