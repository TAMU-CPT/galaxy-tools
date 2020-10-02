#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import sys
import argparse
import copy
import logging
import re
from Bio import SeqIO
logging.basicConfig(level=logging.INFO)


def validate_gbk(genbank_file=None):
    records = SeqIO.parse(genbank_file, "genbank")

    warningCount = 0
    errorCount = 0

    outStr = ""

    for record in records:
        orderPairs = {}
        seenLocus = {}
        pairedRBS = {}
        for feature in record.features:
            if feature.type == "source":
              continue
            if 'locus_tag' not in feature.qualifiers:
              if feature.type != 'regulatory': # Terminators allowed to overlap other features and not have locus_tag, basically anything goes (according to spec provided for this tool)
                outStr += "Warning: Unhandled top-level feature type in " + record.id + " at start location " + str(feature.location.start + 1) + ", type is " + feature.type + "\n"
                warningCount += 1
              continue
                
           # if feature.type == 'CDS' and feature.qualifiers['locus_tag'][0] == "Meda_016":
           #   print((feature.location.parts[0].start))
           #   exit()

            #if feature.type == "mat_peptide":
            #  continue

            if feature.type == "regulatory" and feature.qualifiers['regulatory_class'][0] == 'ribosome_binding_site':
              theKey = str(feature.location.start) + ", " + str(feature.location.end)
              if theKey not in pairedRBS.keys():
                pairedRBS[theKey] = []
              pairedRBS[theKey].append(feature.qualifiers['locus_tag'][0])

            if feature.qualifiers['locus_tag'][0] not in seenLocus.keys():
              seenLocus[feature.qualifiers['locus_tag'][0]] = []
            seenLocus[feature.qualifiers['locus_tag'][0]].append(feature)
        
        for locus in seenLocus.keys():
          geneStart = -1
          geneEnd = -1
          rbsStart = -1
          rbsEnd = -1
          cdsStart = -1
          cdsEnd = -1
          trnaStart = -1
          trnaEnd = -1
          foundIntron = False
          foundSubfeatures = []
          absFeatStart = -1
          absFeatEnd = -1
          for feature in seenLocus[locus]:
            foundSubfeatures.append(feature.type)
            if absFeatStart == -1:
              absFeatStart = feature.location.start
              absFeatEnd = feature.location.end
            else: 
              absFeatStart = min(absFeatStart, feature.location.start)
              absFeatEnd = max(absFeatEnd, feature.location.end)
            if feature.type == 'CDS':
              if cdsStart != -1:
                outStr += "Error: More than one CDS in " + record.id + " at locus_tag " + locus + ".\n"
                errorCount += 1
                cdsStart = min(feature.location.start, cdsStart)
                cdsEnd = max(feature.location.end, cdsEnd)
              else:
                cdsStart = feature.location.start
                cdsEnd = feature.location.end
            elif feature.type == 'regulatory':
              if feature.qualifiers['regulatory_class'][0] == 'ribosome_binding_site':
                if rbsStart != -1:
                  outStr += "Error: More than one RBS in " + record.id + " at locus_tag " + locus + ".\n"
                  errorCount += 1
                  rbsStart = min(feature.location.start, rbsStart)
                  rbsEnd = max(feature.location.end, rbsEnd)
                else:
                  rbsStart = feature.location.start
                  rbsEnd = feature.location.end
              else:
                outStr += "Warning: Unhandled regulator class type in " + record.id + " at locus_tag " + locus + ", type is " + feature.qualifiers['regulatory_class'][0] + ".\n"
                warningCount += 1
            elif feature.type == 'RBS':
              if rbsStart != -1:
                  outStr += "Error: More than one RBS in " + record.id + " at locus_tag " + locus + ".\n"
                  errorCount += 1
                  rbsStart = min(feature.location.start, rbsStart)
                  rbsEnd = max(feature.location.end, rbsEnd)
              else:
                  rbsStart = feature.location.start
                  rbsEnd = feature.location.end
            elif feature.type == "intron":
               startMatch = False
               endMatch = False
               foundIntron = True
               #print(str(feature.location.parts[0].start) + ", " + str(feature.location.parts[0].end))
               for findCDS in seenLocus[locus]:
                 if findCDS.type != "CDS":
                   continue
                 for i in findCDS.location.parts:
                   #print (str(i.start) + ", " + str(i.end))
                   startMatch = startMatch or (i.start == feature.location.parts[0].end)
                   endMatch = endMatch or (i.end == feature.location.parts[0].start)
               if not(startMatch) or not(endMatch):
                 outStr += "Error: Intron feature at locus tag " + locus + " is not bounded by CDS pair.\n"
                 errorCount += 1
            elif feature.type == 'tRNA':
              if trnaStart != -1:
                outStr += "Error: More than one tRNA in " + record.id + " at locus_tag " + locus + ".\n"
                errorCount += 1
                trnaStart = min(feature.location.start, trnaStart)
                trnaEnd = max(feature.location.end, trnaEnd)
              else:
                trnaStart = feature.location.start
                trnaEnd = feature.location.end
            elif feature.type == 'gene':
              if geneStart != -1:
                outStr += "Error: More than one gene with the locus tag " + locus + " in " + record.id + ".\n"
                errorCount += 1
                continue
              else:
                if not(re.search(r'_gt(\d+)', locus)):
                  orderPairs[int(re.findall(r'(\d+)', locus)[-1])] = feature
                geneStart = feature.location.start
                geneEnd = feature.location.end
                
            elif feature.type == "mat_peptide":
                continue
            else:
              outStr += "Warning: Unhandled sub-feature type in " + record.id + " at locus tag " + locus + ", type is " + str(feature.type) + '.\n'
              warningCount += 1

          if geneStart == -1:
              outStr += "Errors (" + str(len(foundSubfeatures)) + " at this location): Subfeatures found for locus_tag " + locus + " but there is no gene feature with an equivalent locus_tag.\n"
              outStr += "\tFeatures " + str(foundSubfeatures) + " at location [" + str(absFeatStart) + ", " + str(absFeatEnd) + "]\n"
              errorCount += len(foundSubfeatures)
              
          elif rbsStart == -1 and cdsStart == -1:
            if trnaStart == -1:
              outStr += "Warning: Gene in " + record.id + " at locus_tag " + locus + ' has no RBS or CDS.\n'
              warningCount += 1
          elif cdsStart == -1:
              outStr += "Warning: Gene in " + record.id + " at locus_tag " + locus + ' has no CDS'
                if foundSubfeatures:
                  outStr += " (Does have the following subfeatures: " + str(foundSubfeatures) + ").\n"
                else:
                  outStr += ".\n"
              warningCount += 1
          elif cdsStart != geneStart and cdsEnd != geneEnd and cdsStart != -1:
              outStr += "Error: CDS in " + record.id + " at locus_tag " + locus + ' does not line up with parent gene.\n'
              outStr += "\tGene: [" + str(geneStart + 1) + ", " + str(geneEnd) + "] --- CDS: [" + str(cdsStart + 1) + ", " + str(cdsEnd) + "]\n"
              errorCount += 1
          #elif rbsStart == -1:
          elif rbsStart != geneStart and rbsEnd != geneEnd and rbsStart != -1:
              shifted = False
              for x in pairedRBS.keys():
                if locus in pairedRBS[x]:
                  shifted = (len(pairedRBS[x]) > 1)
              if not(shifted):
                outStr += "Error: RBS in " + record.id + " at locus_tag " + locus + ' does not line up with parent gene and did not have a possible frameshift equivalent RBS.\n'
                outStr += "\tGene: [" + str(geneStart + 1) + ", " + str(geneEnd) + "] --- RBS: [" + str(rbsStart + 1) + ", " + str(rbsEnd) + "]\n"
                errorCount += 1
          elif ((rbsStart > cdsStart and rbsStart < cdsEnd) or (cdsStart > rbsStart and cdsStart < rbsEnd)) and (rbsStart != -1 and cdsStart != -1):
              outStr += "Error: CDS and RBS overlap in " + record.id + " at locus_tag " + locus + '.\n'
              outStr += "\tCDS: [" + str(cdsStart + 1) + ", " + str(cdsEnd) + "] --- RBS: [" + str(rbsStart + 1) + ", " + str(rbsEnd) + "]\n"
              errorCount += 1
        for x in orderPairs.keys():
          for y in orderPairs.keys():
            if x == y:
              continue
            if x > y and (orderPairs[x].location.start < orderPairs[y].location.start or (orderPairs[x].location.end < orderPairs[y].location.end and orderPairs[x].location.start == orderPairs[y].location.start)):
              outStr += "Error: Gene at locus_tag " + orderPairs[x].qualifiers['locus_tag'][0] + " comes before lower-numbered gene at locus_tag " + orderPairs[y].qualifiers['locus_tag'][0] + "\n"
              outStr += "\t" + orderPairs[x].qualifiers['locus_tag'][0] + ": [" + str(orderPairs[x].location.start + 1) + ", " + str(orderPairs[x].location.end) + "] --- " + orderPairs[y].qualifiers['locus_tag'][0] + ": [" + str(orderPairs[y].location.start + 1) + ", " + str(orderPairs[y].location.end) + "]\n"
              errorCount += 1
        outStr = "QC for " + record.id + " finished with " + str(errorCount) + " errors and " + str(warningCount) + " warnings.\n\n" + outStr
    print(outStr)

if __name__ == "__main__":
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description="Merge GFF3 data into a Genbank file")
    parser.add_argument(
        "genbank_file", type=argparse.FileType("r"), help="Genbank file"
    )
    
    args = parser.parse_args()
    validate_gbk(args.genbank_file)
