#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type, feature_test_true
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    parser.add_argument("changeTo", type=str, help="Feature types to rebase")
    parser.add_argument("changeList", type=str, help="Feature types to create heirarchy with")
    args = parser.parse_args()

    args.changeTo = args.changeTo.split(" ")
    args.changeList = args.changeList.split(" ")

    for record in GFF.parse(args.gff3):
        newRecs = []
        for feature in feature_lambda(record.features, feature_test_true, {}):
            if feature.type in args.changeList:
                if "Parent" in feature.qualifiers.keys():
                    feature.qualifiers["Parent"] = []
                newChain = [feature]
                for x in range(0, len(args.changeTo)):
                    tempFeat = SeqFeature(location=feature.location)
                    tempFeat.type = args.changeTo[len(args.changeTo) - 1 - x]
                    tempFeat.id = feature.id + "_p" + str(len(args.changeTo) - x)
                    tempFeat.ref_db = feature.ref_db
                    tempFeat.ref = feature.ref
                    #tempFeat.sub_features.append(newChain[x])
                    newChain[x].qualifiers["Parent"] = [tempFeat.id]
                    tempFeat.qualifiers["ID"] = [tempFeat.id]
                    tempFeat.qualifiers["Name"] = [feature.qualifiers["Name"][0] + "_p" + str(len(args.changeTo) - x)]
                    newChain.append(tempFeat)
                    
                for x in range(1, len(newChain)):
                    newRecs.append(newChain[x])
                    
        for x in newRecs:
            record.features.append(x)
        GFF.write([record], sys.stdout)