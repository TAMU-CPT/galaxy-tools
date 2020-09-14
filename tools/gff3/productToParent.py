#!/usr/bin/env python
import argparse
import logging
import sys
from cpt_gffParser import gffParse, gffWrite
from gff3 import feature_lambda, feature_test_true, fetchParent
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def extract(qualifier, gff3):
    for record in gffParse(gff3):
        for feature in feature_lambda(record.features, feature_test_true, {}):
            if qualifier in feature.qualifiers:
                tempProd = feature.qualifiers[qualifier][0]
                feature.qualifiers['Name'][0] = feature.qualifiers['Name'][0] + ".p01"
                tempFeat = feature._parent
                tempFeat = tempFeat._parent
                tempFeat.qualifiers[qualifier] = tempProd
                #print("Finding feature:")
                #print(tempFeat)
                #print("-------------------------")
                
        gffWrite([record], sys.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract specified qualifers from features in GFF3', epilog="")
    parser.add_argument('gff3', help='GFF3 File')
    parser.add_argument('qualifier', help='Sepcific qualifier to extract')
    args = parser.parse_args()
    extract(**vars(args))
