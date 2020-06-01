#!/usr/bin/env python

import sys
import argparse
import os
import re
from Bio import SeqIO


class CheckSequence:
    """ 
    SAR endolysin Verification class, which starts with complete FA file, and is shrunk by each function to reveal best candidates of SAR endolysin proteins 
    """


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


    def check_hydrophobicity_and_charge(self,tmd_min=15,tmd_max=20):
        """ verifies the existence of a hydrophobic region within the sequence """
        hydrophobic_residues = "['FIWLVMYCATGS']" # fed through regex
        hits = self.store
        pos_res = "RK"
        neg_res = "DE"

        if self.size > 50:
            seq = self.seq[0:50]
        else:
            seq = self.seq 
        for tmd_size in range(tmd_min, tmd_max, 1):
            for i in range(0,len(seq)-tmd_size,1):
                tmd_seq = str(seq[i:i+tmd_size]) 
                if re.search((hydrophobic_residues+"{"+str(tmd_size)+"}"),tmd_seq):
                    charge_seq = str(seq[:i]) # look at all of the residues up to the index of the sequence that we're currently at
                    charge = charge_check(charge_seq,pos_res,neg_res)
                    storage_dict(self=self,tmd_size=tmd_size,tmd_seq=tmd_seq,hits=hits,charge_seq=charge_seq,charge=charge,perc_cont=percent_calc(tmd_seq,"GA",float(tmd_size)))
                elif "K" in tmd_seq[0] and re.search((hydrophobic_residues+"{"+str(tmd_size-1)+"}"),tmd_seq[1:]): # check frontend snorkels
                    charge_seq = str(seq[:i]) # look at all of the residues up to the index of the sequence that we're currently at
                    charge = charge_check(charge_seq,pos_res,neg_res)
                    storage_dict(self=self,tmd_size=tmd_size,tmd_seq=tmd_seq,hits=hits,charge_seq=charge_seq,charge=charge,perc_cont=percent_calc(tmd_seq,"GA",float(tmd_size)))
                elif "K" in tmd_seq[-1] and re.search((hydrophobic_residues+"{"+str(tmd_size-1)+"}"),tmd_seq[:-1]): # check backend snorkels
                    charge_seq = str(seq[:i]) # look at all of the residues up to the index of the sequence that we're currently at
                    charge = charge_check(charge_seq,pos_res,neg_res)
                    storage_dict(self=self,tmd_size=tmd_size,tmd_seq=tmd_seq,hits=hits,charge_seq=charge_seq,charge=charge,perc_cont=percent_calc(tmd_seq,"GA",float(tmd_size)))
                continue
        
        return hits


### Extra "helper" functions
def storage_dict(self,tmd_size,tmd_seq,hits,charge_seq,charge,perc_cont): # probably not good to call "self" a param here...definitley not PEP approved...
    """ organize dictionary for hydrophobicity check """
    if self.name not in hits:
        hits[self.name] = {}
        hits[self.name]["description"] = str(self.description)
        hits[self.name]["sequence"] = str(self.seq)
        hits[self.name]["size"] = str(self.size)
        #GAcont = str((str(self.seq).count("G")+str(self.seq).count("A"))/int(self.size)*100)
        #hits[self.name]["GAcont"] = "{:.2f}%".format(float(GAcont))
        if "TMD_"+str(tmd_size) not in hits[self.name]:
            hits[self.name]["TMD_"+str(tmd_size)] = []
            hits[self.name]["TMD_"+str(tmd_size)].append([tmd_seq,charge_seq,charge,perc_cont])
        else:
            hits[self.name]["TMD_"+str(tmd_size)].append([tmd_seq,charge_seq,charge,perc_cont])
    else:
        if "TMD_"+str(tmd_size) not in hits[self.name]:
            hits[self.name]["TMD_"+str(tmd_size)] = []
            hits[self.name]["TMD_"+str(tmd_size)].append([tmd_seq,charge_seq,charge,perc_cont])
        else:
            hits[self.name]["TMD_"+str(tmd_size)].append([tmd_seq,charge_seq,charge,perc_cont]) 


def percent_calc(sequence,residues,size):
    """ Calculate the percent of a set of residues within an input sequence """
    counted = {}
    for aa in sequence:
        #print(aa)
        if aa in counted:
            counted[aa] += 1
        else:
            counted[aa] = 1
    residue_amt = 0
    for res_of_interest in residues:
        try:
            residue_amt += counted[res_of_interest]
        except KeyError:
            residue_amt += 0
    ratio = residue_amt/size
    
    return ratio*100


def charge_check(charge_seq,pos_res,neg_res):
    charge = 0
    for aa in charge_seq:
        if aa in pos_res:
            charge += 1
        if aa in neg_res:
            charge -= 1
    return charge

if __name__ == "__main__":
    pass
