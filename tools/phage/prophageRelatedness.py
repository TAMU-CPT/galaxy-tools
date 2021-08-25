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


def parseXML(blastxml, outFile): # Modified from intron_detection
    blast = []
    for iter_num, blast_record in enumerate(NCBIXML.parse(blastxml), 1):
        align_num = 0
        outFile.write("Query ID\tQuery Length\tTotal Number of Hits\n")
        outFile.write("%s\t%d\t%d\n\n" % (blast_record.query_id, blast_record.query_length, len(blast_record.alignments)))
        
        for alignment in blast_record.alignments:

            align_num += 1
            gi_nos = str(alignment.accession)
            blast_gene = []
            for hsp in alignment.hsps:

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
                        "align_len": hsp.align_length,
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
        
def compPhage(inRec, outFile, padding = 1.2, cutoff = .3, numReturn = 20):
    inRec = parseXML(inRec, outFile)
    res = []
    for group in inRec:
      window = floor(padding * float(group[0]["query_length"]))
      group = sorted(group, key = lambda x: x["sbjct_range"][0])
      hspGroups = []
      lastInd = len(res)
      
      for x in range(0, len(group)):
        hspGroups.append([group[x]])
        startBound = group[x]["sbjct_range"][0]
        endBound = startBound + window
        for hsp in group[x + 1:]:
          if hsp["sbjct_range"][0] >= startBound and hsp["sbjct_range"][1] <= endBound:
            hspGroups[-1].append(hsp)
        
      for x in disjointSets(superSets(hspGroups)):
        res.append(x)
      
      maxID = 0.0
      for x in res[lastInd:]:
        sumID = 0.0
        totAlign = 0
        for y in x:
          totAlign += y["align_len"] 
          sumID += float(y["identity"])
        x.append(totAlign)
        x.append(sumID / float(x[0]["query_length"]))
        maxID = max(maxID, x[-1])
        
    res = sorted(res, key = lambda x: x[-1], reverse = True)

    outList = []
    outNum = 0
    for x in res:
        if outNum == numReturn or x[-1] < cutoff:
          break
        outNum += 1
        outList.append(x)
        
          

        
#    Original request was that low scoring clusters would make it to the final results IF
#    they were part of an Accession cluster that did have at least one high scoring member.

    
    outFile.write("Accession Number\tCluster Start Location\tEnd Location\tSubject Cluster Length\t# HSPs in Cluster\tTotal Aligned Length\t% of Query Aligned\tOverall % Query Identity\tOverall % Subject Identity\tComplete Accession Info\n")
    for x in outList:
      minStart = min(x[0]["sbjct_range"][0], x[0]["sbjct_range"][1])
      maxEnd = max(x[0]["sbjct_range"][0], x[0]["sbjct_range"][1])
      if "|gb|" in x[0]["match_id"]:
        startSlice = x[0]["match_id"].index("gb|") + 3
        endSlice = (x[0]["match_id"][startSlice:]).index("|")
        accOut = x[0]["match_id"][startSlice: startSlice + endSlice]
      else:
        accOut = x[0]["gi_nos"]
      for y in x[0:-2]:
# ("\t%.3f\t" % (x[-1]))
        minStart = min(minStart, y["sbjct_range"][0])
        maxEnd = max(maxEnd, y["sbjct_range"][1])
      outFile.write(accOut + "\t" + str(minStart) + "\t" + str(maxEnd) + "\t" + str(maxEnd - minStart + 1) + "\t" + str(len(x) - 1) + "\t" + str(x[-2]) + ("\t%.3f" % (float(x[-2]) / float(x[0]['query_length']) * 100.00)) + ("\t%.3f" % (x[-1] * 100.00)) + ("\t%.3f" % (float(x[-2]) / float(maxEnd - minStart + 1) * 100.00)) + "\t" + x[0]["match_id"] + "\n")
   
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
        "--cutoff",
        help="Gap minimum (Default -1, set to a negative number to allow overlap)",
        default=.3,
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
