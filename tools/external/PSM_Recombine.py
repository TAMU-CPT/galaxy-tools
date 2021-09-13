#!/usr/bin/env python
import argparse
import logging
from Bio import SeqIO

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify shine-dalgarno sequences")
    parser.add_argument("psmTable", type=argparse.FileType("r"))
    parser.add_argument('gbkList', type=argparse.FileType("r"), nargs="+")
    args = parser.parse_args()
    
    gbkRecs = []
    recIDs = []
    recFlatten = [] # Can only seek argparse file once

    for f in args.gbkList:
      tempRecs = SeqIO.parse(f, "genbank")
      for rec in tempRecs:
        recFlatten.append(rec)

    for line in args.psmTable:
      lineElems = line.split("\t")
      numGenes = 0
      accession = ""
      lineOut = ""
      if recIDs == []:
        for i in lineElems:
          recIDs.append(i.strip())
          lineOut += i.strip() + "\t"
          for rec in recFlatten:
              if i.strip() in rec.id or rec.id in i.strip():
                gbkRecs.append(rec)
        lineOut += "No. of phages in which gene is present\tBest Database Match"
        print(lineOut)
        continue
     
      for i in range(0, len(lineElems)):
        checkFeat = lineElems[i].strip()
        if checkFeat == "-":
          lineOut += "(-)\t"
          continue
        else:
          lineOut += checkFeat + "\t"
        numGenes += 1
        if accession == "":
          for feat in gbkRecs[i].features:
            if "locus_tag" in feat.qualifiers.keys() and feat.qualifiers["locus_tag"][0] == checkFeat:
              if "protein_id" in feat.qualifiers.keys():
                accession = feat.qualifiers["protein_id"][0]
                break # Comment out if we need to get more info
      lineOut += str(numGenes) + "\t" + accession
      print(lineOut)


