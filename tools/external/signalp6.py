#!/usr/bin/env python

import os
import argparse
import subprocess
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Python shell to suppress pipe errors in SignalP6")

    parser.add_argument("--organism", type=str)
    parser.add_argument("--inFile", type=str)

    args = parser.parse_args()
    cmd = "signalp6 -org " + args.organism + " -od ./subDir --format png --mode slow-sequential -ff " + args.inFile
    try :
        subprocess.run(cmd, stderr=sys.stdout.buffer, shell=True) # Just to keep any sub-module errors Prophet throws to keep Galaxy from thinking the job failed
    except :                            # The except is irrelevant, just needs to complete
        res = 0
    
    
    

