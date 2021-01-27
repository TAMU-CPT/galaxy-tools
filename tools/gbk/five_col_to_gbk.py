#!/usr/bin/env python
import argparse
import logging
import sys
from Bio import SeqIO, Alphabet
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

def fiveColToGbk(tabIn, seqIn):
  seqList = list(SeqIO.parse(seqIn, "fasta"))
  seqOut = None
  startLoc = 0
  endLoc = -1
  featList = []
  featType = ""
  recOut = []
  for line in tabIn:
    print(line)
    if line[0] == ">":
      if featType:
        #if len(seqList) == 1
        #  seqOut = seqList[0].seq
        #else:
        for x in seqList:
          if x.id == recID:
            seqOut = x.seq   
            seqOut.alphabet = IUPAC.IUPACUnambiguousDNA()         
        if not seqOut:
          log.error("Unable to find associated sequence for 5-column record " + recID + ", unable to construct Genbank")   
        recOut.append(SeqRecord(seqOut, id=recID, features=sorted(featList, key=lambda x: x.location.start)))
        featType = ""
        featList = []
        seqOut = None
      recID = line[9:].strip()
      continue
  
    fields = line.split("\t")
    if len(fields) > 2 and fields[2]:
      if featType:
        featList.append(SeqFeature(featLoc, featType, qualifiers=featQuals))
      featType = fields[2][0:-1]
      featLoc = FeatureLocation(int(fields[0]), int(fields[1]), strand=0)
      featQuals = {}
    elif fields[0]:
      featLoc.parts.append(FeatureLocation(int(fields[0]), int(fields[1]), strand=0))
    elif fields[3]:
      if fields[3] in featQuals.keys():
        featQuals[fields[3]].append(fields[4][0:-1])
      else:
        featQuals[fields[3]] = [fields[4][0:-1]]
    # else blank line

  if len(seqList) == 1:
    seqOut = seqList[0].seq
    seqOut.alphabet = IUPAC.IUPACUnambiguousDNA() 
  else:
    for x in seqList:
      if x.id == recID:
        seqOut = x.seq
        seqOut.alphabet = IUPAC.IUPACUnambiguousDNA() 

  if not seqOut:
    log.error("Unable to find associated sequence for 5-column record " + recID + ", unable to construct Genbank")
   
  recOut.append(SeqRecord(seqOut, id=recID, features=sorted(featList, key=lambda x: x.location.start)))
  
  return recOut    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a Genbank file into five column format"
    )
    parser.add_argument("tabIn", type=argparse.FileType("r"), help="Five Column tabular input")
    parser.add_argument("seqIn", type=argparse.FileType("r"), help="Associated Sequence(s) (Fasta input)")

    args = vars(parser.parse_args())
    for rec in fiveColToGbk(**args):
      SeqIO.write(rec, sys.stdout, "genbank")
