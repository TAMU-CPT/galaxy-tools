#!/usr/bin/env python
import sys
import copy
import argparse
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def lipoP_gff(lipoIn, gff3In):

    orgIDs = {}
    orgID = ""

    # Take and parse the txt output into a sequence of records
    # Dict of X records, with the ID as key and an array Y of each cleavage site as the value, 
    for row in lipoIn:
       if row.startswith('#'):
           orgID = ""
           continue
 
       rowElem = row.split('\t')

       orgID = rowElem[0]
       if not (orgID in orgIDs.keys()): 
           orgIDs[orgID] = []

       if rowElem[2] == "CleavII": 
           orgIDs[orgID].append(int(rowElem[3]))#, int(rowElem[4])))

    # Rebase
    for gff in GFF.parse(gff3In):
        for xRec in gff.features: 
            findCleave = ""
            for cleID in orgIDs.keys():
                if cleID == xRec.id:
                    findCleave = cleID
                    break
            if findCleave == "":
                continue

            i = 0
            for cleaveBase in orgIDs[findCleave]:
                tempQuals = xRec.qualifiers.copy()
                i += 1
                tempQuals['ID'] = xRec.id + "_cleavage_" + str(i)
                xRec.sub_features.append(SeqFeature(FeatureLocation(xRec.location.start + cleaveBase, xRec.location.start + cleaveBase + 1), type="cleavage_site",  strand = xRec.location.strand, qualifiers = tempQuals))
            
        GFF.write([gff], sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='add parent gene features to CDSs')
    parser.add_argument('lipoIn', type=argparse.FileType("r"), help='LipoP tool\'s .txt output')
    parser.add_argument('gff3In', type=argparse.FileType("r"), help='GFF3 to rebase LipoP results')
    args = parser.parse_args()
    lipoP_gff(**vars(args))
