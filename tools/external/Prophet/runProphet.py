#!/usr/bin/env python

import os
import argparse
import subprocess


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CPT's very own modified Efetch")

    parser.add_argument("fasIn", type=str) # current place holder until I determine how best to use the current user's email from Galaxy

    parser.add_argument("gffIn", type=str)

    args = parser.parse_args()
    cmd = "/galaxy/tools/cpt2/galaxy-tools/tools/external/Prophet/ProphET_standalone.pl --fasta" + args.fasIn + " --gff_in" + args.gffIn + "--outdir ./tempOut"
    try :
        subprocess.run(cmd, shell=True)#, stdout=output)
    except :
        res = 0
    
    
    return

