#!/usr/bin/env python
import sys
import copy
import argparse
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from gff3 import feature_lambda, feature_test_type, get_id

def lipoP_gff(lipoIn, gff3In, jBrowseOut):

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

       if rowElem[2] == "CleavII":
           if not (orgID in orgIDs.keys()): 
               orgIDs[orgID] = [] 
           orgIDs[orgID].append(int(rowElem[3]))#, int(rowElem[4])))

    # Rebase
    for gff in GFF.parse(gff3In):
        keepSeq = []
        for xRec in gff.features:
            cdss = list(feature_lambda(xRec.sub_features, feature_test_type, {'type': 'CDS'}, subfeatures=False))
            findCleave = ""
            cdsOff = 0
            for cds in cdss:
                if cds.id in orgIDs:
                    findCleave = cds.id
                    break
                cdsOff += 1
            if findCleave == "":
                if not jBrowseOut:
                    keepSeq.append(xRec)
                continue

            if jBrowseOut:
                xRec.sub_features = []

            i = 0
            for cleaveBase in orgIDs[findCleave]:
                tempQuals = xRec.qualifiers.copy()
                i += 1
                tempQuals['ID'] = xRec.id + "_cleavage_" + str(i)
                
                xRec.sub_features.append(SeqFeature(FeatureLocation(cdss[cdsOff].location.start + (cleaveBase * 3) - 1, cdss[cdsOff].location.start + (cleaveBase * 3) + 1), type="cleavage_site",  strand = xRec.location.strand, qualifiers = tempQuals))
            keepSeq.append(xRec)

        gff.features = keepSeq    
        GFF.write([gff], sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='add parent gene features to CDSs')
    parser.add_argument('lipoIn', type=argparse.FileType("r"), help='LipoP tool\'s .txt output')
    parser.add_argument('gff3In', type=argparse.FileType("r"), help='GFF3 to rebase LipoP results')
    parser.add_argument('--jBrowseOut', type=bool, default = False, help='Prepare Output for jBrowse instance')
    args = parser.parse_args()
    lipoP_gff(**vars(args))
