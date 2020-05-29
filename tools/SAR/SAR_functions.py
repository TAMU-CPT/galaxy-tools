#!/usr/bin/env python

import sys
import argparse
import os
import re
from Bio import SeqIO


class CheckSequence:
    """ SAR endolysin Verification class, which starts with complete FA file, and is shrunk by each function to reveal best candidates of SAR endolysin proteins """

    def __init__(self, protein_name, protein_data):
        self.name = protein_name
        self.seq = protein_data.seq
        self.description = protein_data.description
        self.size = len(self.seq)
        self.store = {}

    def check_sizes(self,min,max):
        """ check the minimum and maximum peptide lengths """
        if self.size < min:
            print("too small")
        elif self.size > max:
            print("too large")
        else:
            print(f"{self.name} : {self.seq}")
            return True

    def check_hydrophobicity(self,tmd_min=15,tmd_max=20):
        """ verifies the existence of a hydrophobic region within the sequence """
        hydrophobic_residues = "['FIWLVMYCATGS']" # fed through regex
        hits = self.store
        if self.size > 50:
            seq = self.seq[0:50]
        else:
            seq = self.seq 
        for tmd_size in range(tmd_min, tmd_max, 1):
            for i in range(0,len(seq)-tmd_size,1):
                check_seq = str(seq[i:i+tmd_size]) # fix this tomorrow
                if re.search((hydrophobic_residues+"{"+str(tmd_size)+"}"),check_seq):
                    storage_dict(self=self,tmd_size=tmd_size,check_seq=check_seq,hits=hits)
                elif "K" in check_seq[0] and re.search((hydrophobic_residues+"{"+str(tmd_size-1)+"}"),check_seq[1:]): # check frontend snorkels
                    storage_dict(self=self,tmd_size=tmd_size,check_seq=check_seq,hits=hits)
                elif "K" in check_seq[-1] and re.search((hydrophobic_residues+"{"+str(tmd_size-1)+"}"),check_seq[:-1]): # check backend snorkels
                    storage_dict(self=self,tmd_size=tmd_size,check_seq=check_seq,hits=hits)
                continue
        
        return hits


### Extra "helper" functions
def storage_dict(self,tmd_size,check_seq,hits):
    if self.name not in hits:
        hits[self.name] = {}
        hits[self.name]["description"] = str(self.description)
        hits[self.name]["sequence"] = str(self.seq)
        hits[self.name]["size"] = str(self.size)
        if "TMD_"+str(tmd_size) not in hits[self.name]:
            hits[self.name]["TMD_"+str(tmd_size)] = []
            hits[self.name]["TMD_"+str(tmd_size)].append([check_seq])
        else:
            hits[self.name]["TMD_"+str(tmd_size)].append([check_seq])
    else:
        if "TMD_"+str(tmd_size) not in hits[self.name]:
            hits[self.name]["TMD_"+str(tmd_size)] = []
            hits[self.name]["TMD_"+str(tmd_size)].append([check_seq])
        else:
            hits[self.name]["TMD_"+str(tmd_size)].append([check_seq]) 

if __name__ == "__main__":
    pass
