#!/usr/bin/env python
import sys
import argparse
from BCBio import GFF
from Bio.Blast import NCBIXML

class IntronFinder(object):
    pass
    # parse xml into appropriate data structure
    # def __init__():

    # checker function(s)
    # def check_strand():
    # def check_distance():

    # look through gi nos, find matching blast hits
    # call checker functions to make sure intron is legit
    # delete records from data structure that are not introns
    # def find_introns():

    # merge 2 or more seq records into one
    # def modify_gff3(?):

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Intron detection')
    parser.add_argument('genes', type=file, help='GFF3 Gene Calls')
    parser.add_argument('blastp', type=file, help='blast XML protein results')
    args = parser.parse_args()
