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
    startLoc='0',
    endLoc='1'
):

    
    if locusMode:
      for record in SeqIO.parse(genbank_file, "genbank"):
        record.features = sorted(record.features, key=lambda x: x.location.start)
        startPos = 0
        endPos = 1
        if endLoc == "":
          endLoc = record.features[-1].qualifiers['locus_tag'][0] # Attempt to get last feature if none supplied
        featOut = []
        addFeats = False
        nextEnd = False
        for feature in record.features:
          if 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] == startLoc and not addFeats:
            addFeats = True
            startPos = feature.location.start
          elif addFeats and 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] == endLoc:
            nextEnd = True
          elif nextEnd and (('locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] != endLoc) or 'locus_tag' not in feature.qualifiers):  # inclusive end (To include final feature if desired)
            addFeats = False  # Attempts to get entire feature tree of last feature

          if addFeats:
            if endPos < feature.location.end:
              endPos = feature.location.end
            featOut.append(feature)
        

        if revCom:
          finSeq = (record.seq[startPos: endPos]).reverse_complement()
        else:
          finSeq = record.seq[startPos: endPos] 
        for x in featOut:
          x.location = FeatureLocation(x.location.start - startPos, x.location.end - startPos, x.location.strand)
        yield [
                    SeqRecord(
                        Seq(str(finSeq).strip(), record.seq.alphabet),
                        id=record.id,
                        features=featOut,
                    )
              ]
    else:
      for record in SeqIO.parse(genbank_file, "genbank"):
        featOut = []
        for feature in record.features:
          if feature.location.start >= int(startLoc) - 1 and feature.location.end < int(endLoc):
                featOut.append(feature)
        if revCom:
          finSeq = (record.seq[int(startLoc) - 1: int(endLoc)]).reverse_complement()
        else:
          finSeq = record.seq[int(startLoc) - 1: int(endLoc)]
        for x in featOut:
          x.location = FeatureLocation(x.location.start - (int(startLoc) - 1), x.location.end - (int(startLoc) - 1), x.location.strand)
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
