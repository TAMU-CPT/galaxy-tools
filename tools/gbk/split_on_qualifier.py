#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def splitByQual(gbkIn=None, key=[], val=[], baseAhead=0, baseBehind=0):
    res = []
    
    for record in SeqIO.parse(gbkIn, "genbank"):
        contStart = -1
        contEnd = -1
        prevStrand = 1
        outFeat = []

        for x in sorted(record.features, key=lambda x: x.location.start):
            for i in key:
              if i in list(x.qualifiers.keys()) and any(j in val for j in x.qualifiers[i]):
                if contStart < 0:
                  contStart = x.location.start
                  contEnd = x.location.end
                  prevStrand = x.location.strand  
                  outFeat.append(x)
                elif (x.location.start >= contStart and x.location.start <= contEnd) or (x.location.end >= contStart and x.location.end <= contEnd):
                  contStart = min(x.location.start, contStart)
                  contEnd = max(x.location.end, contEnd)
                  if prevStrand != x.location.strand:
                    prevStrand = 0
                  outFeat.append(x)
                else:
                  if prevStrand < 0:
                    contEnd = min(contEnd + baseAhead, len(record.seq) - 1)
                    contStart = max(contStart - baseBehind, 0)
                  else:
                    contEnd = min(contEnd + baseBehind, len(record.seq) - 1)
                    contStart = max(contStart - baseAhead, 0)
                  outSeq = record.seq[contStart:contEnd]
                  res.append(SeqRecord(outSeq, record.id, record.name, record.description, record.dbxrefs, sorted(outFeat, key=lambda x: x.location.start), record.annotations, record.letter_annotations))
                  contStart = x.location.start
                  contEnd = x.location.end
                  prevStrand = x.location.strand  
                  outFeat = [x]                  
        if contStart != -1:
                  if prevStrand < 0:
                    contEnd = min(contEnd + baseAhead, len(record.seq) - 1)
                    contStart = max(contStart - baseBehind, 0)
                  else:
                    contEnd = min(contEnd + baseBehind, len(record.seq) - 1)
                    contStart = max(contStart - baseAhead, 0)
                  outSeq = record.seq[contStart:contEnd]
                  res.append(SeqRecord(outSeq, record.id, record.name, record.description, record.dbxrefs, sorted(outFeat, key=lambda x: x.location.start), record.annotations, record.letter_annotations))
                  contStart = -1
                  contEnd = -1
                  prevStrand = 1

    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset a genbank file")
    parser.add_argument('gbkIn', type=argparse.FileType("r"), help="Multi-record Genbank file")
    parser.add_argument('--key', type=str, default="locus_tag", required=False, help="Qualifier to split on")
    parser.add_argument('--val', type=str, default="locus_tag", required=False, help="Value to match qualifier against")
    parser.add_argument("--baseAhead", type=int, help="Number of extra bases upstream to extract", default=0)
    parser.add_argument("--baseBehind", type=int, help="Number of extra bases downstream to extract", default=0)

    args = parser.parse_args()
    args.key = args.key.split(" ")
    args.val = args.val.split(" ")

    for record in splitByQual(**vars(args)):
        SeqIO.write(record, sys.stdout, "genbank")
