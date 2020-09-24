#!/usr/bin/env python
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Section off a subset of a Genbank file", epilog=""
    )
    parser.add_argument(
        "inFile", type=argparse.FileType("r"), help="Tabular file"
    )
    parser.add_argument(
        "gffFile", type=argparse.FileType("w"), help="Genbank file"
    )
    parser.add_argument(
        "fastaFile", type=argparse.FileType("w"), help="Tabular file"
    )
    args = vars(parser.parse_args())

    writeMode = 1
    outGFF = ""
    outFASTA = ""
    for line in args['inFile'].readlines():
      if line == "###\n":
        continue
      elif line == "##FASTA\n":
        writeMode = 2
      elif writeMode == 1:
        args['gffFile'].write(line)
      elif writeMode == 2:
        args['fastaFile'].write(line)
      
#    args['gffFile'].write(outGFF)
#    args['fastaFile'].write(outFASTA)
