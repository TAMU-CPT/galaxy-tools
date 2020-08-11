#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import sys
import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

def makeSubset(
    genbank_file=None,
    locusMode=False,
    revCom=False,
    startLoc='',
    endLoc=''
):

    numStart = -1
    numEnd = -1
    lastEnd = 0
    
    if locusMode:
      for record in SeqIO.parse(genbank_file, "genbank"):
        record.features = sorted(record.features, key=lambda x: x.location.start)
        for feature in record.features:
          lastEnd = int(max(lastEnd, max(feature.location.start, feature.location.end)))
          if 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] == startLoc:
            if numStart == -1:
              numStart = int(min(feature.location.start, feature.location.end))
            else:
              numStart = int(min(numStart, min(feature.location.start, feature.location.end)))

          elif 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] == endLoc:
            numEnd = int(max(numEnd, max(feature.location.start, feature.location.end)))
      if startLoc == '':
        numStart = 0
      if endLoc == '':
        numEnd = lastEnd
      
      if numStart == -1:
        exit(2)
      if numEnd == -1:
        exit(2)
    else:  
      if startLoc == '':
        numStart = 0
      else:
        numStart = int(startLoc) - 1
      if endLoc == '':
        numEnd = lastEnd
      else:
        numEnd = int(endLoc) - 1
          

    
    for record in SeqIO.parse(genbank_file, "genbank"):
        featOut = []
        for feature in record.features:
          if feature.location.start >= numStart and feature.location.end < numEnd:
                featOut.append(feature)
        if revCom:
          finSeq = (record.seq[numStart: numEnd]).reverse_complement()
        else:
          finSeq = record.seq[numStart: numEnd]
        for x in featOut:
          if revCom:
            x.location = FeatureLocation(numEnd - (x.location.end - numStart), numEnd - (x.location.start - numStart), x.location.strand)
          else:
            x.location = FeatureLocation(x.location.start - numStart, x.location.end - numStart, x.location.strand)
        featOut = sorted(featOut, key=lambda x: x.location.start)
        yield [
                    SeqRecord(
                        Seq(str(finSeq).strip(), record.seq.alphabet),
                        id=record.id,
                        features=featOut,
                    )
              ]


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Section off a subset of a Genbank file", epilog=""
    )
    parser.add_argument(
        "genbank_file", type=argparse.FileType("r"), help="Genbank file"
    )
    parser.add_argument(
        "--locusMode", action="store_true", help="Use locus tags"
    )
    parser.add_argument("--revCom", action="store_true", help="Reverse complement sequence")
    parser.add_argument("--startLoc", help="Starting Location", default='0')
    parser.add_argument("--endLoc", help="Ending Location", default='1')
    args = vars(parser.parse_args())
    for seq in makeSubset(**args):
        SeqIO.write(seq, sys.stdout, "genbank")
