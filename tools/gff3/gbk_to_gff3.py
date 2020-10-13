#!/usr/bin/env python

import argparse
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from cpt_gffParser import gffSeqFeature, gffWrite

def makeGffFeat(inFeat, num, recName):
    if "codon_start" in inFeat.qualifiers.keys():
      shift = int(inFeat.qualifiers["codon_start"][0]) - 1
    else:
      shift = "."
    if "locus_tag" in inFeat.qualifiers.keys():
      name = inFeat.qualifiers["locus_tag"][0] + "." + inFeat.type 
      if num > 0:
        name += "." + str(num)
    else:
      name = recName + "." + inFeat.type + "." + str(num)
    
    outFeat = gffSeqFeature(inFeat.location, inFeat.type, '', inFeat.strand, name, inFeat.qualifiers, None, None, None, shift, 0, "GbkToGff")
    outFeat.qualifiers["ID"] = [name]  
    return outFeat

def main():
    parser = argparse.ArgumentParser( description='Convert GenBank flat files to GFF3 format')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GBK file' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output GFF file to be created' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    features_skipped_count = 0
    outRec = []

    failed = 0

    for rec in SeqIO.parse(open(args.input_file, "r"), "genbank"):
        recID = rec.name

        if len(str(rec.seq)) > 0:
            seqs_pending_writes = True
            outSeq = str(rec.seq)
            seqLen = len(outSeq)

        locBucket = {}
        outFeats = []
        topTypeDict = {}
        for feat in rec.features:
            if "locus_tag" not in feat.qualifiers.keys(): #Allow metadata features and other features with no ID (Output warning?) - AJC
              if feat.type not in topTypeDict.keys():
                topTypeDict[feat.type] = 1
              else:
                topTypeDict[feat.type] += 1
              outFeats.append(makeGffFeat(feat, topTypeDict[feat.type], recID))
              continue
            elif feat.qualifiers["locus_tag"][0] not in locBucket.keys():
              locBucket[feat.qualifiers["locus_tag"][0]] = []
            locBucket[feat.qualifiers["locus_tag"][0]].append(feat)

        for locus in locBucket.keys():
          if len(locBucket[locus]) == 1: # No heirarchy to be made
            outFeats.append(makeGffFeat(locBucket[locus][0], 0, recID))
            continue
          topFeat = None
          midFeat = None
          bottomFeats = []
          typeDict = {}
          minLoc = locBucket[locus][0].location.start
          maxLoc = locBucket[locus][0].location.end
          for feat in locBucket[locus]:
            # If we want to make our own top-level feat?
            minLoc = min(minLoc, feat.location.start)
            maxLoc = max(maxLoc, feat.location.end)
            
            # Gene->mRNA->CDS included as example, to add other feature-heirarchys in the appropriate slot
            if feat.type in ['gene']:
              if not topFeat:
                topFeat = feat
              # Else handle multiple top features
            elif feat.type in ['mRNA', 'tRNA', 'rRNA']:
              if not midFeat:
                midFeat = feat
              # Else handle multiple mid feats (May need another elif type-in-list statement if we actually expect a list of mid feats)
            else:
              if feat.type not in typeDict.keys():
                typeDict[feat.type] = 1
              else:
                typeDict[feat.type] += 1 
              bottomFeats.append(feat)
            
          #if not topFeat: # Make our own top-level feature based off minLoc, maxLoc bounds
            
          for x in typeDict.keys(): # If only 1, set it to 0 so we don't append a number to the name
            if typeDict[x] == 1:    # Else, set to 1 so that we count up as we encounter the features
              typeDict[x] = 0
            else:
              typeDict[x] = 1
          
          if not topFeat:
            sys.stderr.write("Unable to create a feature heirarchy at location [%d, %d] with features: \n" % (minLoc, maxLoc))
            for x in locBucket[locus]:
              sys.stderr.write(str(x))
              sys.stderr.write('\n')
              failed = 1
            continue

          outFeats.append(makeGffFeat(topFeat, 0, recID))
          if midFeat: # Again, need a new if statement if we want to handle multiple mid-tier features
              outFeats[-1].sub_features.append(makeGffFeat(midFeat, 0, recID))
              outFeats[-1].sub_features[-1].qualifiers["Parent"] = [outFeats[-1].id]
              for x in bottomFeats:
                outFeats[-1].sub_features[-1].sub_features.append(makeGffFeat(x, typeDict[x.type], recID))
                outFeats[-1].sub_features[-1].sub_features[-1].qualifiers["Parent"] = [outFeats[-1].sub_features[-1].id]
                typeDict[x.type] += 1
          else: # No midFeat, append bottom feats directly to top feats (Make our own midFeat?)
              for x in bottomFeats:
                outFeats[-1].sub_features.append(makeGffFeat(x, typeDict[x.type], recID))
                outFeats[-1].sub_features[-1].qualifiers["Parent"] = [outFeats[-1].id]
                typeDict[x.type] += 1
        outRec.append(SeqRecord(rec.seq, recID, rec.name, rec.description, rec.dbxrefs, outFeats, rec.annotations, rec.letter_annotations)) 

    gffWrite(outRec, ofh)    
    exit(failed) # 0 if all features handled, 1 if unable to handle some


if __name__ == '__main__':
    main()








