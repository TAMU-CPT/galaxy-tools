#!/usr/bin/env python

from cpt_gffParser import gffParse, gffWrite
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO, SeqFeature, SeqRecord
from Bio.Seq import Seq
from Bio.Data import CodonTable
from gff3 import feature_lambda, feature_test_true
import csv
import argparse
from cpt_gffParser import gffParse, gffWrite

def codonCheck(thisFeat):
  featOut = thisFeat
  if thisFeat.type != "CDS":
    featOut.shift = '.'
  elif "codon_start" in thisFeat.qualifiers.keys():
    featOut.shift = int(thisFeat.qualifiers["codon_start"][0]) - 1
  else:
    featOut.shift = 0  
  for x in range(0, len(featOut.sub_features)):
    resFeat = codonCheck(featOut.sub_features[x])
    featOut.sub_features[x] = resFeat
  
  return featOut 
  
def translationCheck(thisFeat, coSeq=None, table = 11):
  featOut = thisFeat
  finDict = {}
  for x in thisFeat.qualifiers.keys():
    if x == "translation":
      if coSeq != None:
        finTranslate = ""
        for loc in thisFeat.location.parts:
          if featOut.strand == -1:
            shiftSeq = (coSeq[loc.start: loc.end - featOut.shift]).reverse_complement()
          else:
            shiftSeq = coSeq[loc.start + featOut.shift: loc.end]
          try:
            finTranslate += str(shiftSeq.translate(table=table, cds=True))
          except CodonTable.TranslationError as cte:
          #  log.info("Translation issue at %s: %s", record.id, cte)
            finTranslate += str(shiftSeq.translate(table=table, cds=False))
        
        finDict["translation"] = [finTranslate]
    else:
      finDict[x] = thisFeat.qualifiers[x]
  featOut.qualifiers = finDict
  for x in range(0, len(featOut.sub_features)):
    resFeat = translationCheck(featOut.sub_features[x], coSeq, table)
    featOut.sub_features[x] = resFeat
  
  return featOut 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fixes poorly constructed GFF from Genbank to GFF tool"
    )
    parser.add_argument("gff3In", type=argparse.FileType("r"), help="GFF3 source file")
    parser.add_argument(
        "--fasta",
        type=argparse.FileType("r"),
        help="Fasta File for translating"
    )
    parser.add_argument(
        "--table",
        type=int,
        default=0,
        help="Translation table to use",
        choices=range(0, 23),
    )
    args = parser.parse_args()
    
    if args.table > 0:
      seq_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
    else:
      seq_dict = {}
    outRec = []
    for record in list(gffParse(args.gff3In, seq_dict)):
      for ind in range(0, len(record.features)):
        tempFeat = codonCheck(record.features[ind])
        qualFeat = tempFeat
        if args.table > 0:
          tempFeat = translationCheck(qualFeat, record.seq, args.table)
        else:
          tempFeat = translationCheck(qualFeat)
        record.features[ind] = tempFeat
      outRec.append(record)
    gffWrite(outRec)
