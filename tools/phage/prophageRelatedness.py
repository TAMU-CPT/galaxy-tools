#!/usr/bin/env python
import sys
import re
import itertools
import argparse
import hashlib
import copy
from math import floor
from Bio.Blast import NCBIXML
from collections import OrderedDict
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()


def parseXML(blastxml): # Modified from intron_detection
    blast = []
    for iter_num, blast_record in enumerate(NCBIXML.parse(blastxml), 1):
        align_num = 0
        print(dir(blast_record))
        for alignment in blast_record.alignments:
            print(dir(alignment))
            align_num += 1
            gi_nos = str(alignment.accession)
            blast_gene = []
            for hsp in alignment.hsps:
                print(dir(hsp))
                x = float(hsp.identities - 1) / ((hsp.query_end) - hsp.query_start)
                nice_name = blast_record.query

                if " " in nice_name:
                    nice_name = nice_name[0 : nice_name.index(" ")]
                blast_gene.append(
                    {
                        "gi_nos": gi_nos,
                        "sbjct_length": alignment.length,
                        "query_length": blast_record.query_length,
                        "sbjct_range": (hsp.sbjct_start, hsp.sbjct_end),
                        "query_range": (hsp.query_start, hsp.query_end),
                        "name": nice_name,
                        "evalue": hsp.expect,
                        "identity": hsp.identities,
                        "identity_percent": x,
                        "hit_num": align_num,
                        "iter_num": iter_num,
                        "match_id": alignment.title.partition(">")[0],
                    }
                )
            blast.append(blast_gene)
            
        
    return blast


def test_true(feature, **kwargs):
    return True

def superSets(inSets):
    inSets.sort(key=len, reverse = True)
    nextInd = 0
    res = []
    for i in range(0, len(inSets)):
      if i == 0:
        res.append(inSets[i])
        continue
      for par in res:
        complete = True
        for x in inSets[i]:
          if not (x in par): 
            complete = False
        if complete:
          break # Subset of at least one member
      if not complete:
        res.append(inSets[i])
    return res

def disjointSets(inSets):
    inSets.sort(key = lambda x: x[0]["sbjct_range"][0])
    res = [inSets[0]]
    for i in range(1, len(inSets)):
       disjoint = True
       for elem in inSets[i]:
         for cand in res:
           if elem in cand:
             disjoint = False
             break
         if not disjoint:
           break
       if disjoint:
         res.append(inSets[i])
    return res
        
def compPhage(inRec, outFile, padding = 1.2, numReturn = 20):
    inRec = parseXML(inRec)
    res = []
    for group in inRec:
      window = floor(padding * float(group[0]["query_length"]))
      group = sorted(group, key = lambda x: x["sbjct_range"][0])
      hspGroups = []
      for x in range(0, len(group)):
        hspGroups.append([group[x]])
        startBound = group[x]["sbjct_range"][0]
        endBound = startBound + window
        for hsp in group[x + 1:]:
          if hsp["sbjct_range"][0] >= startBound and hsp["sbjct_range"][1] <= endBound:
            hspGroups[-1].append(hsp)
        
      res.append(disjointSets(superSets(hspGroups)))
      
      maxID = 0.0
      for x in res[-1]:
        sumID = 0.0
        for y in x:
          sumID += float(y["identity"])
        x.append(sumID / float(x[0]["query_length"]))
        maxID = max(maxID, x[-1])
      res[-1].append(maxID)
        
    res = sorted(res, key = lambda x: x[-1], reverse = True)

    outList = []
    outNum = 0
    for x in res:
      for y in x[0:-1]:
        outNum += 1
        outList.append(y)
        if outNum == numReturn:
          break
      if outNum == numReturn:
        break

    outList.sort(key = lambda x: x[-1], reverse = True)

    outFile.write("Accession Number\tScore\tCluster Start Location\tEnd Location\t# HSPs in Cluster\tComplete Accession Info\n")

    for x in outList:
      minStart = x[0]["sbjct_range"][0]
      maxEnd = x[0]["sbjct_range"][1]
      for y in x[0:-1]:
        minStart = min(minStart, y["sbjct_range"][0])
        maxEnd = max(maxEnd, y["sbjct_range"][1])
      outFile.write(x[0]["gi_nos"] + "\t" + str(x[-1]) + "\t" + str(minStart) + "\t" + str(maxEnd) + "\t" + str(len(x) - 1) + "\t" + x[0]["match_id"] + "\n")
   
    #accession start end number

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Intron detection")
    parser.add_argument(
        "inRec", type=argparse.FileType("r"), help="blast XML protein results"
    )
    parser.add_argument(
        "--outFile",
        type=argparse.FileType("w"),
        help="Output Error Log",
        default="./compPro.tsv",
    )
    parser.add_argument(
        "--padding",
        help="Gap minimum (Default -1, set to a negative number to allow overlap)",
        default=1.2,
        type=float,
    )
    parser.add_argument(
        "--numReturn",
        help="Gap maximum in genome (Default 10000)",
        default=20,
        type=int,
    )
    args = parser.parse_args()

    compPhage(**vars(args))
