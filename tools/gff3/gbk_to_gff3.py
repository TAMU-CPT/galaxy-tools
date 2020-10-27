#!/usr/bin/env python

import argparse
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from cpt_gffParser import gffSeqFeature, gffWrite

bottomFeatTypes = ["exon", "RBS", "CDS"]

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
        seekingParent = []
        for feat in rec.features:
            if "locus_tag" not in feat.qualifiers.keys(): #Allow metadata features and other features with no ID (Output warning?) - AJC
              if feat.type in bottomFeatTypes:
                seekingParent.append([feat, [], []]) # [Feature, all parent candidates, strongest parent candidates]
                continue
              elif feat.type not in topTypeDict.keys():
                topTypeDict[feat.type] = 1
              else:
                topTypeDict[feat.type] += 1
              outFeats.append(makeGffFeat(feat, topTypeDict[feat.type], recID))
              continue
            elif feat.qualifiers["locus_tag"][0] not in locBucket.keys():
              locBucket[feat.qualifiers["locus_tag"][0]] = []
            locBucket[feat.qualifiers["locus_tag"][0]].append(feat)

        for locus in locBucket.keys():
          minLoc = locBucket[locus][0].location.start
          maxLoc = locBucket[locus][0].location.end
          for feat in locBucket[locus]:
            minLoc = min(minLoc, feat.location.start)
            maxLoc = max(maxLoc, feat.location.end)
          for x in seekingParent:
            if x[0].location.start >= minLoc and x[0].location.end <= maxLoc:
              x[1].append(locus)
            if x[0].location.start == minLoc or x[0].location.end == maxLoc:
              x[2].append(locus)

        for x in seekingParent: #Reformat to [Feature, Locus, Unused/Free]
          if len(x[2]) == 1:
            finList = ""
            if len(x[1]) > 1:
              for loc in x[1]:
                if loc != x[2][0]:
                  finList += loc + ", "
              finList = str(x[0].type) + " had no locus tag set in .gbk file, automatically derived. Other, weaker candidate(s) were " + finList[0:-2] + "."
            else:
              finList = str(x[0].type) + " had no locus tag set in .gbk file, automatically derived."
            if "Notes" not in x[0].qualifiers.keys():
              x[0].qualifiers["Notes"] = []
            x[0].qualifiers["Notes"].append(finList)
            x[1] = x[2][0]
          elif len(x[2]) > 1:
            candidate = x[2][0] #Arbitrarily choose first one
            finList = ""
            strongList = ""
            for loc in x[2]:
              if loc != candidate:
                finList += loc + ", "
                strongList += loc + ", "
            for loc in x[1]:
              if loc not in x[2]:
                finList += loc + ", "
            finList = str(x[0].type) + " had no locus tag set in .gbk file, automatically derived. Other candidate(s) were " + finList[0:-2] + " (Equally strong candidate(s): " + strongList[0:-2] + ")."
            if "Notes" not in x[0].qualifiers.keys():
              x[0].qualifiers["Notes"] = []
            x[0].qualifiers["Notes"].append(finList)
            x[1] = candidate
          elif len(x[1]) == 1:
            x[1] = x[1][0]
            if "Notes" not in x[0].qualifiers.keys():
              x[0].qualifiers["Notes"] = []
            finList = str(x[0].type) + " had no locus tag set in .gbk file, automatically derived."
            x[0].qualifiers["Notes"].append(finList)
          elif len(x[1]) > 1:
            candidate = x[1][0] #Arbitrarily choose first one
            finList = ""
            for loc in x[1]:
              if loc != candidate:
                finList += loc + ", "
            finList = str(x[0].type) + " had no locus tag set in .gbk file, automatically derived. Other candidates were " + finList[0:-2] + "."
            if "Notes" not in x[0].qualifiers.keys():
              x[0].qualifiers["Notes"] = []
            x[0].qualifiers["Notes"].append(finList)
            x[1] = candidate
          else:
            sys.stderr.write("Warning: Unable to find potential parent for feature with no locus tag of type " + str(x[0].type) + " at location [" + str(x[0].location.start) + ", " + str(x[0].location.end) + "].\n")
            if x[0].type not in topTypeDict.keys():
              topTypeDict[x[0].type] = 1
            else:
              topTypeDict[x[0].type] += 1
            outFeats.append(makeGffFeat(x[0], topTypeDict[x[0].type], recID))
            
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
          
          for x in seekingParent:
            if type(x[1]) != "list" and locus == x[1]:
              x[0].qualifiers["locus_tag"] = [locus]
              bottomFeats.append(x[0])
              if x[0].type not in typeDict.keys():
                typeDict[x[0].type] = 1
              else:
                typeDict[x[0].type] += 1 
              
          
            
              
            
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








