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


    def check_hydrophobicity_and_charge(self,sar_min=15,sar_max=20,perc_residues="SGA"):
        """ verifies the existence of a hydrophobic region within the sequence """
        hydrophobic_residues = "['FIWLVMYCATGS']" # fed through regex
        hits = self.store
        pos_res = "RK"
        neg_res = "DE"

        if self.size > 50:
            seq = self.seq[0:50]
        else:
            seq = self.seq 
        for sar_size in range(sar_min, sar_max, 1):
            for i in range(0,len(seq)-sar_size,1):
                sar_seq = str(seq[i:i+sar_size])
                if re.search((hydrophobic_residues+"{"+str(sar_size)+"}"),sar_seq):
                    charge_seq, charge, perc_cont, sar_coords, nterm_coords, cterm_coords, sar_start = rep_funcs(self,seq,i,pos_res,neg_res,sar_seq,perc_residues,sar_size)
                    storage_dict(self=self,sar_size=sar_size,sar_seq=sar_seq,hits=hits,charge_seq=charge_seq,charge=charge,perc_cont=perc_cont,nterm_coords=nterm_coords,sar_coords=sar_coords,cterm_coords=cterm_coords,sar_start=sar_start)
                    #print("TMDSIZE: {}\tINDEX: {}".format(sar_size,i+1))
                elif "K" in sar_seq[0] and re.search((hydrophobic_residues+"{"+str(sar_size-1)+"}"),sar_seq[1:]): # check frontend snorkels
                    charge_seq, charge, perc_cont, sar_coords, nterm_coords, cterm_coords, sar_start = rep_funcs(self,seq,i,pos_res,neg_res,sar_seq,perc_residues,sar_size)
                    storage_dict(self=self,sar_size=sar_size,sar_seq=sar_seq,hits=hits,charge_seq=charge_seq,charge=charge,perc_cont=perc_cont,nterm_coords=nterm_coords,sar_coords=sar_coords,cterm_coords=cterm_coords,sar_start=sar_start)
                    #print("TMDSIZE: {}\tINDEX: {}".format(sar_size,i+1))
                elif "K" in sar_seq[-1] and re.search((hydrophobic_residues+"{"+str(sar_size-1)+"}"),sar_seq[:-1]): # check backend snorkels
                    charge_seq, charge, perc_cont, sar_coords, nterm_coords, cterm_coords, sar_start = rep_funcs(self,seq,i,pos_res,neg_res,sar_seq,perc_residues,sar_size)
                    storage_dict(self=self,sar_size=sar_size,sar_seq=sar_seq,hits=hits,charge_seq=charge_seq,charge=charge,perc_cont=perc_cont,nterm_coords=nterm_coords,sar_coords=sar_coords,cterm_coords=cterm_coords,sar_start=sar_start)
                    #print("TMDSIZE: {}\tINDEX: {}".format(sar_size,i+1))
                continue
        
        return hits

    def shrink_results(self,sar_min=15,sar_max=20,perc_residues="SGA"):
        """ removes repetiive hits, keeps only the shortest and longest of each SAR domain """
        compare_candidates = {}
        hits = self.check_hydrophobicity_and_charge()
        for sar_name, data in hits.items():
            print(sar_name)
            compare_candidates[sar_name] = {}
            #print("\nThese are the values: {}".format(v))
            #count_of_times = 0
            tmd_log = []
            for sar_size in range(sar_max,sar_min-1,-1):
                if "TMD_"+str(sar_size) in data:
                    tmd_log.append(sar_size)
                    #print(tmd_log) 
                    #print(data["TMD_"+str(sar_size)])
                    for idx,the_data in enumerate(data["TMD_"+str(sar_size)]):
                        #print(f"This is the index: {idx}")
                        #print(f"This is the list of data at this index: {the_data}")
                        if the_data[-1] in compare_candidates[sar_name]:
                            compare_candidates[sar_name][the_data[-1]]["count"] += 1
                            compare_candidates[sar_name][the_data[-1]]["size"].append(sar_size)
                            compare_candidates[sar_name][the_data[-1]]["index"].append(idx)
                        else:
                            compare_candidates[sar_name][the_data[-1]] = {}
                            compare_candidates[sar_name][the_data[-1]]["count"] = 1 
                            compare_candidates[sar_name][the_data[-1]]["size"] = [sar_size]
                            compare_candidates[sar_name][the_data[-1]]["index"] = [idx]
            hits[sar_name]["biggest_sar"] = tmd_log[0]
            print("Biggest sar --> "+str(hits[sar_name]["biggest_sar"]))
        for sar_name, compare_data in compare_candidates.items():
            for data in compare_data.values():
                if len(data["size"]) >= 3:
                    #print(f"{each_size} --> {data}")
                    minmax = [min(data["size"]),max(data["size"])]
                    nonminmax = [x for x in data["size"] if x not in minmax]
                    nonminmax_index = []
                    for each_nonminmax in nonminmax:
                        v = data["size"].index(each_nonminmax)
                        x = data["index"][v]
                        nonminmax_index.append(x)
                    nons = zip(nonminmax,nonminmax_index)
                    for value in nons:
                        #hits[sar_name]["TMD_"+str(value[0])] = hits[sar_name]["TMD_"+str(value[0])].pop(value[1])
                        hits[sar_name]["TMD_"+str(value[0])][value[1]] = [""]

        return hits


def rep_funcs(self,seq,loc,pos_res,neg_res,sar_seq,perc_residues,sar_size):
    """ run a set of functions together before sending the results to the storage dictionary """

    charge_seq = str(seq[:loc])
    charge = charge_check(charge_seq,pos_res,neg_res)
    perc_cont = percent_calc(sar_seq,perc_residues,int(sar_size))
    sar_start = loc
    sar_coords = "{}..{}".format(loc,loc+sar_size)
    nterm_coords = "{}..{}".format("0",loc-1)
    cterm_coords = "{}..{}".format(loc+sar_size+1,self.size)

    return charge_seq, charge, perc_cont, sar_coords, nterm_coords, cterm_coords, sar_start


### Extra "helper" functions
def storage_dict(self,sar_size,sar_seq,hits,charge_seq,charge,perc_cont,nterm_coords,sar_coords,cterm_coords,sar_start): # probably not good to call "self" a param here...definitley not PEP approved...
    """ organize dictionary for hydrophobicity check """
    if self.name not in hits:
        hits[self.name] = {}
        hits[self.name]["description"] = str(self.description)
        hits[self.name]["sequence"] = str(self.seq)
        hits[self.name]["size"] = str(self.size)
        #GAcont = str((str(self.seq).count("G")+str(self.seq).count("A"))/int(self.size)*100)
        #hits[self.name]["GAcont"] = "{:.2f}%".format(float(GAcont))
        if "TMD_"+str(sar_size) not in hits[self.name]:
            hits[self.name]["TMD_"+str(sar_size)] = []
            hits[self.name]["TMD_"+str(sar_size)].append([sar_seq,charge_seq,charge,perc_cont,nterm_coords,sar_coords,cterm_coords,sar_start])
        else:
            hits[self.name]["TMD_"+str(sar_size)].append([sar_seq,charge_seq,charge,perc_cont,nterm_coords,sar_coords,cterm_coords,sar_start])
    else:
        if "TMD_"+str(sar_size) not in hits[self.name]:
            hits[self.name]["TMD_"+str(sar_size)] = []
            hits[self.name]["TMD_"+str(sar_size)].append([sar_seq,charge_seq,charge,perc_cont,nterm_coords,sar_coords,cterm_coords,sar_start])
        else:
            hits[self.name]["TMD_"+str(sar_size)].append([sar_seq,charge_seq,charge,perc_cont,nterm_coords,sar_coords,cterm_coords,sar_start]) 


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
    my_ratios = []
    for res_of_interest in residues:
        try:
            residue_amt = counted[res_of_interest]
        except KeyError:
            residue_amt = 0
        ratio = residue_amt/size
        my_ratios.append((round(ratio*100,2)))
    
    res_rat = zip(residues,my_ratios)

    return res_rat


def charge_check(charge_seq,pos_res,neg_res):
    charge = 0
    for aa in charge_seq:
        if aa in pos_res:
            charge += 1
        if aa in neg_res:
            charge -= 1
    return charge

if __name__ == "__main__":
    sequence = "MAGBYYYTRLCVRKLRKGGGHP"
    residues = "YL"
    size = len(sequence)
    print(size)
    v = percent_calc(sequence,residues,size)
    print(v)
    for i in v:
        print(i)

