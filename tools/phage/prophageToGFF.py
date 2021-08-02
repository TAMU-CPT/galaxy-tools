#!/usr/bin/env python
# vim: set fileencoding=utf-8
import os
import argparse
from Bio import SeqIO
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="pat")

def convertTable(prophageIn, fastaIn):
    disallowArray = ["&", ",", ";", "="]
    validArray = ["%26", "%2C", "%3B", "%3D"]
    print("##gff-version 3")
    seqSource = SeqIO.parse(fastaIn, "fasta")
    for rec in seqSource:
        orgName = rec.id
    for row in prophageIn:
        line = row.strip().split("\t") 
        if len(line) != 9 or row.startswith("Accession Number"):
            continue
        
        geneName = line[0]
        geneStart = min(int(line[1]), int(line[2]))
        geneEnd = max(int(line[1]), int(line[2]))
        geneQuals ="ID=%s;Name=%s;Score=%s;Notes=" % (geneName, geneName, line[7])
        for x in line[8]:
          try:
            insIndex = disallowArray.index(x)
            geneQuals = geneQuals + validArray[insIndex]
          except:
            geneQuals = geneQuals + x
        geneQuals = geneQuals + validArray[1] + " contains %s HSPs;" % line[4]
        print("%s\tRelatedProphages\tgene\t%s\t%s\t.\t+\t.\t%s" % (orgName, geneStart, geneEnd, geneQuals))
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="rebase gff3 features against parent locations", epilog=""
    )
    parser.add_argument(
        "prophageIn", type=argparse.FileType("r"), help="Parent GFF3 annotations"
    )
    parser.add_argument("fastaIn", type=argparse.FileType("r"), help="Genome Sequence")

    args = parser.parse_args()

    convertTable(**vars(args))
