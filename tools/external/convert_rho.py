from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio import SeqIO, SeqFeature, SeqRecord
from Bio.Seq import Seq, UnknownSeq
from cpt_gffParser import gffSeqFeature, gffWrite
import sys

gffOut = open(sys.argv[1], "w")
res = []
for recPath in sys.argv[2:]:
  recName = recPath[25:-4]
  termList = []
  termNum = 0
  maxLoc = 0
  recFile = open(recPath, "r")
  for line in recFile:
    line = line.strip()
    if line[0:23] == "PREDICTED REGION NUMBER":
      if termNum != 0:
        termList.append(gffSeqFeature(FeatureLocation(startLoc, endLoc, strand), "terminator", "", strand, termID, {"Note": notes, "ID": [termID]}, source = "RhoTermPredict"))
      termNum += 1
      termID = recName + "_T" + str(termNum)
      if line[-9:-1] == "POSITIVE":
        strand = 1
      else:
        strand = -1
      continue
    elif line[0:16] == "Genomic sequence":
      indDash = line.index("-")
      indComma = line.index(",")
      indClose = line.index(" ):")
      startLoc = int(line[53:indDash]) - 1
      endLoc = int(line[indDash+1:indComma])
      maxLoc = max(maxLoc, endLoc)
      indCG = line.index("c/g")
      notes = [line[indCG:indClose]]
      continue
    elif line[0:21] == "Palindromic sequences" or line[0:15] == "PAUSE-CONSENSUS":
      notes.append(line.strip())
      continue
  termList.append(gffSeqFeature(FeatureLocation(startLoc, endLoc, strand), "terminator", "", strand, termID, {"Note": notes, "ID": [termID]}, source = "RhoTermPredict"))
  res.append(SeqRecord.SeqRecord(UnknownSeq(maxLoc), id=recName, name=recName, description="RhoTermPredict on " + recName, features=termList))
gffWrite(res, gffOut)

    
