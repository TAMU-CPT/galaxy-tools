#!/usr/bin/env python

from cpt_gffParser import gffParse, gffWrite
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO, SeqFeature, SeqRecord
from Bio.Seq import Seq
from gff3 import feature_lambda, feature_test_true
import csv
import argparse
from cpt_gffParser import gffParse, gffWrite



def checkParentageBoundaries(thisFeat):
  errMess = ""
  errCount = 0
  warnMess = ""
  warnCount = 0
  if str(thisFeat.type).upper() == "GENE":
    foundMRNA = False
    for x in thisFeat.sub_features:
      if str(x.type).upper() == "MRNA":
        foundMRNA = True
        if x.location.start != thisFeat.location.start or x.location.end != thisFeat.location.end:
          errMess += "Error: MRNA " + x.id + " does not match location of parent gene " + thisFeat.id + " (Expected both start and end to match for both features).\n"
          errCount += 1
      elif str(x.type).upper() != "REGULATORY" and str(x.type).upper() != "SEQUENCE_SECONDARY_STRUCTURE":
        errMess += "Error: Gene " + thisFeat.id + " has unexpected subfeature of type " + x.type + " (Expected subtypes for gene: mRNA, regulatory, sequence_secondary_structure).\n"
        errCount += 1
    if not foundMRNA:
      warnMess += "Warning: No MRNA subfeature found for gene " + thisFeat.id + ".\n" 
      warnCount += 1
  elif str(thisFeat.type).upper() == "MRNA":
    hasCDS = False
    hasExon = False
    cdsBound = False
    exonBoundStart = False
    exonBoundEnd = False
    exon1 = None
    exon2 = None
    for x in thisFeat.sub_features:
      if str(x.type).upper() == "CDS":
        hasCDS = True
        if x.location.start < thisFeat.location.start or x.location.end > thisFeat.location.end:
          errMess += "Error: CDS " + x.id + " falls outside range of parent feature " + thisFeat.id + " (mRNA: [" + str(thisFeat.location.start) + ", " + str(thisFeat.location.end) + "], CDS [" + str(x.location.start) + ", " + str(x.location.end) + "]).\n"
          errCount += 1
        if x.location.start == thisFeat.location.start or x.location.end == thisFeat.location.end:
          cdsBound = True
      elif str(x.type).upper() == "EXON":
        hasExon = True
        if not exon1:
          exon1 = x.location
        elif not exon2:
          exon2 = x.location
        if x.location.start < thisFeat.location.start or x.location.end > thisFeat.location.end:
          errMess += "Error: Exon " + x.id + " falls outside range of parent feature " + thisFeat.id + " (mRNA: [" + str(thisFeat.location.start) + ", " + str(thisFeat.location.end) + "], Exon [" + str(x.location.start) + ", " + str(x.location.end) + "]).\n"
          errCount += 1
        if x.location.start == thisFeat.location.start:
          exonBoundStart = True
        if x.location.end == thisFeat.location.end:
          exonBoundEnd = True
      elif str(x.type).upper() != "INTRON" and str(x.type).upper() != "SHINE_DALGARNO_SEQUENCE":
        errMess += "Error: mRNA " + thisFeat.id + " has unexpected subfeature of type " + x.type + " (Expected subtypes for mRNA: CDS, exon, intron, and Shine_Dalgarno_sequence).\n"
        errCount += 1

    if hasCDS:
      if hasExon and not (exonBoundStart and exonBoundEnd):
        if not cdsBound:
          errMess += "Error: mRNA " + thisFeat.id + " has no combination of CDSs or Exons that cover both the start and end.\n"
          errCount += 1
      elif not cdsBound:
        errMess += "Error: mRNA " + thisFeat.id + " has CDS subfeatures but none that cover either the start or the end of the feature.\n" 
        errCount += 1
    elif hasExon and not (exonBoundStart and exonBoundEnd):
      errMess += "Error: mRNA " + thisFeat.id + " has no Exon(s) that cover both the start and the end of the feature.\n"
      errCount += 1
    elif not foundCDS and not foundExon:
      warnMess += "Warning: No CDS or Exon subfeatures found for mRNA " + thisFeat.id + ".\n" 
      warnCount += 1
    if exon1 and exon2:
      if not(exon1.start > exon2.end or exon2.start > exon1.end):
        errMess += "Error: Overlapping exons in mRNA " + thisFeat.id + ".\n"
  elif thisFeat.sub_features and (str(thisFeat.type).upper() == "CDS" or str(thisFeat.type).upper() == "EXON" or str(thisFeat.type).upper() == "INTRON" or str(thisFeat.type).upper() == "SHINE_DALGARNO_SEQUENCE"):
    errMess += "Error: " + thisFeat.type + " " + thisFeat.id + " has subfeatures (Expected to be final, bottom-level feature).\n"
    errCount += 1

  for x in thisFeat.sub_features:
    resErr, resErrNum, resWarn, resWarnNum = checkParentageBoundaries(x)
    errMess += resErr
    errCount += resErrNum
    warnMess += resWarn
    warnCount += resWarnNum
  return errMess, errCount, warnMess, warnCount  

  
def tooManyFeatCheck(thisFeat):
  mrnaCount = 0
  exonCount = 0
  errMess = ""
  errCount = 0
  for x in thisFeat.sub_features:
    if str(x.type).upper() == "MRNA":
      mrnaCount += 1
    elif str(x.type).upper() == "EXON":
      exonCount += 1 

  if mrnaCount > 1:
    errMess += "Error: Too many mRNA features in " + str(thisFeat.id) + ", expected 1, found " + str(mrnaCount) + ".\n"
    errCount += 1
  if exonCount > 2:
    errMess += "Error: Too many exon features in " + str(thisFeat.id) + ", expected a maximum of 2, found " + str(exonCount) + ".\n"
    errCount += 1
  
  for x in thisFeat.sub_features:
    resMess, resCount = tooManyFeatCheck(x)
    errMess += resMess
    errCount += resCount
  
  return errMess, errCount

def dupCheck(thisFeat):
  warnMess = ""
  warnCount = 0
  if len(thisFeat.location.parts) > 1:
    warnMess += "Warning: Two features with ID " + str(thisFeat.id) + " found, ensure that join-type location is intended.\n"
    warnCount += 1
  
  for x in thisFeat.sub_features:
    resMess, resCount = dupCheck(x)
    warnMess += resMess
    warnCount += resCount
  
  return warnMess, warnCount

def qualCheck(thisFeat, autoFix=False):
  # The old bad-character checks are now built into the parser, this only checks Genbank conversion problems
  # Genbank qualifiers are delinieated by double qoutes. If a user wishes to use double qoutes within a Note or similar field,
  # the Genbank spec for escaping double qoute is two double qoutes in a row.
  warnOut = ""
  warnCount = 0
  featOut = thisFeat
  for x in thisFeat.qualifiers.keys():
    for sentInd in range(0, len(thisFeat.qualifiers[x])):
      sentence = thisFeat.qualifiers[x][sentInd]
      newSent = ""
      skip = False
      changed = False
      for ind in range(0, len(sentence)):
        if skip:
          continue
        newSent += sentence[ind]
        if sentence[ind] == '"':
          if ind == len(sentence) - 1 or sentence[ind + 1] != '"':
            if autoFix:
              newSent += '"'
              changed = True
              warnOut += "Warning: Qualifier " + x + " in feature " + thisFeat.id + " modified to be GFF-to-Genbank compliant.\n"
              warnCount += 1
            else:
              warnOut += "Warning: Qualifier " + x + " in feature " + thisFeat.id + " has a double-qoute character (\") at position " + str(ind) + ", which will break the GFF-to-Genbank conversion. If conversion to Genbank is desired, place another \" immediately after it to meet the Genbank spec.\n"
              warnCount += 1
          else: # Already compliant
            skip = True 
      if changed:       
        featOut.qualifiers[x][sentInd] = newSent  
  for x in range(0, len(featOut.sub_features)):
    resMess, resCount, resFeat = qualCheck(featOut.sub_features[x], autoFix)
    warnOut += resMess
    warnCount += resCount
    featOut.sub_features[x] = resFeat
  if warnOut and autoFix:
    return warnOut, warnCount, featOut 
  return warnOut, warnCount, thisFeat   
          



def table_annotations(gff3In, out_errorlog, autoFix = True, removeAnnote = []):

    outRec = []
    
    for record in list(gffParse(gff3In)):
        numError = 0
        numWarning = 0
        errorMessage = ""
        warnMessage = ""
        topFeats = []
        metaFeats = []
        for x in record.features:
          if x.type != "annotation" and x.type != "contig":
            topFeats.append(x)
          else:
            metaFeats.append(x)
        topFeats.sort(key=lambda x: x.location.start)

        for featInd in range(0, len(metaFeats)):
          resetID = False
          feat = metaFeats[featInd]
          if feat.id == "" or feat.id == None:
            warnMessage += "Warning: No ID found for " + feat.type + " feature at ["+ str(feat.location.start) + ", " + str(feat.location.end) + "].\n"
            numWarning += 1
            feat.id = "at location ["+ str(feat.location.start) + ", " + str(feat.location.end) + "] (No ID)"
            resetID = True
          resWarn, resNum, featOut = qualCheck(feat, autoFix)
          warnMessage += resWarn
          numWarning += resNum
          if resNum > 0 and autoFix:
            if resetID:
              featOut.id = None
            metaFeats[featInd] = featOut

        for featInd in range(0, len(topFeats)):
          feat = topFeats[featInd]
          if feat.id == "" or feat.id == None:
            warnMessage += "Warning: No ID found for " + feat.type + " feature at ["+ str(feat.location.start) + ", " + str(feat.location.end) + "].\n"
            numWarning += 1
            feat.id = "at location ["+ str(feat.location.start) + ", " + str(feat.location.end) + "] (No ID)"

          resWarn, resNum, featOut = qualCheck(feat, autoFix)
          warnMessage += resWarn
          numWarning += resNum
          if resNum > 0 and autoFix:
            if resetID:
              featOut.id = None
            topFeats[featInd] = featOut

          if str(feat.type).upper() == "CDS" or str(feat.type).upper() == "EXON" or str(feat.type).upper() == "INTRON" or str(feat.type).upper() == "SHINE_DALGARNO_SEQUENCE" or str(feat.type).upper() == "MRNA":
            errorMessage += "Error: Bad top-level feature " + feat.id + " of type " + feat.type + " at ["+ str(feat.location.start) + ", " + str(feat.location.end) + "].\n"
            numError += 1

          resErr, resNum, resWarn, resWarnNum = checkParentageBoundaries(feat)
          errorMessage += resErr
          numError += resNum
          warnMessage += resWarn
          numWarning += resWarnNum
          resErr, resNum = tooManyFeatCheck(feat)
          errorMessage += resErr
          numError += resNum
          resWarn, resNum = dupCheck(feat)
          warnMessage += resWarn
          numWarning += resNum
          
 
          checkInd = featInd + 1
          while checkInd < len(topFeats):
            if feat.location.end < topFeats[checkInd].location.start:
              checkInd = len(topFeats)
            else:
              if feat.type == "tRNA" or topFeats[checkInd].type == "tRNA":
                allowOver = 10
              else:
                allowOver = 60
              badTop = not((feat.type == "gene" or feat.type == "tRNA") and (topFeats[checkInd].type == "gene" or topFeats[checkInd].type == "tRNA"))
              if allowOver < abs(feat.location.end - topFeats[checkInd].location.start):
                errorMessage += "Error: " + feat.type + " " + feat.id +" overlaps " + topFeats[checkInd].type + " " + topFeats[checkInd].id + " by " + str(abs(feat.location.end - topFeats[checkInd].location.start)) + " bases, maximum of " + str(allowOver) + " allowed for " + feat.type + " and " + topFeats[checkInd].type
                if badTop:
                  errorMessage += " (Note: Probably caused by bad top-level feature).\n"
                else:
                  errorMessage += ".\n"
                numError += 1
              checkInd += 1

        ##
        out_errorlog.write("=======================================================================\n")  
        out_errorlog.write("Validation for " + record.id + " finished with " + str(numError) + " errors and " + str(numWarning) + " warnings.\n")
        out_errorlog.write("=======================================================================\n")
        out_errorlog.write(errorMessage)
        out_errorlog.write(warnMessage)

        
    # For each terminator/tRNA
    # Bother Cory later

        if removeAnnote and not autoFix:
          record.annotations = {}
          outFeats = []
          for x in record.features:
            if x.type in removeAnnote:
              outFeats.append(x)
          record.features = outFeats

        if autoFix:
          finFeats = topFeats + metaFeats
          if removeAnnote:
            outFeats = []
            for x in finFeats:
              if x.type not in removeAnnote:
                outFeats.append(x)
            record.features = outFeats
          else:
            record.features = finFeats
          
          
        if autoFix or removeAnnote:
          record.features.sort(key=lambda x: x.location.start)
          outRec.append(record)

    if autoFix or removeAnnote:
      gffWrite(outRec)

#  print(dir(record.features))
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate that contents of GFF3 are correct"
    )
    parser.add_argument("gff3In", type=argparse.FileType("r"), help="GFF3 source file")
    parser.add_argument(
        "--out_errorlog",
        type=argparse.FileType("w"),
        help="Output Error Log",
        default="test-data/errorlog.txt",
    )
    parser.add_argument(
        "--autoFix", action="store_true", help="Fix qualifiers to be genbank compatible"
    )
    parser.add_argument(
        "--removeAnnote", type=str, help="Remove Annotations", default=""
    )
    args = parser.parse_args()
    if args.removeAnnote:
      args.removeAnnote = args.removeAnnote.split(" ")
    
    table_annotations(**vars(args))
