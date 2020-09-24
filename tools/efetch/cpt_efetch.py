#!/usr/bin/env python

import sys
from time import sleep 
import os
from os import path
from Bio import Entrez
from Bio import SeqIO
from urllib.error import HTTPError
import argparse
from helperFunctions import awk_files, is_dir


#Entrez.email = "curtisross@tamu.edu"

class CPTEfetch:
    """ Object that has built in functions to retrieve data from NCBI. Initially constructued to retreive GB and FA files from the nuccore and protein NCBI databases """

    def __init__(self,email,acc,db,ret_type):
        self.email = email
        self.acc = acc
        self.db = db
        self.ret_type = ret_type
        Entrez.email = self.email


    def __repr__(self):
        return "<accession: {} | database: {} | return_type: {}>".format(self.acc,self.db,self.ret_type)


    def retrieve_data(self,sleep_time=5):
        if self.ret_type == "genbank":
            try:
                net_handle = Entrez.efetch(db=self.db,id=self.acc,rettype="gbwithparts",retmode="text")
            except HTTPError:
                sleep(sleep_time)
                net_handle = Entrez.efetch(db=self.db,id=self.acc,rettype="gbwithparts",retmode="text")
        else:
            try:
                net_handle = Entrez.efetch(db=self.db,id=self.acc,rettype="fasta",retmode="text")
            except HTTPError:
                sleep(sleep_time)
                net_handle = Entrez.efetch(db=self.db,id=self.acc,rettype="fasta",retmode="text")
        return net_handle.read()


    def write_record(self,name,st,galaxy=True):
        record = self.retrieve_data(sleep_time=st)
        if galaxy:
            """with name as f:
                print(name)
                f.write(record)"""
            #with open(f"{name}_{str(self.acc)}.{str(self.ret_type)}","w") as file:
            with open(str(name)+"_"+str(self.acc)+"."+str(self.ret_type),"w") as file:
                file.write(record)
        else:
            #with open(f"{str(self.acc)}.{str(self.ret_type)}","w") as file:
            with open(str(self.acc)+"."+str(self.ret_type),"w") as file:
                file.write(record)


if __name__ == "__main__":
        ##### Arguments
    parser = argparse.ArgumentParser(description="CPT's very own modified Efetch")

    parser.add_argument("--email",
                        type=str,
                        help="Entrez Required Email") # current place holder until I determine how best to use the current user's email from Galaxy

    parser.add_argument("--input",
                        type=str,
                        required=True,
                        action="append",
                        #nargs="*",
                        help='accession input"')

    parser.add_argument("--db",
                        type=str,
                        required=True,
                        choices=("protein", "nuccore"),
                        help="choose protein or nuccore database to do query")

    parser.add_argument("--ret_format",
                        type=str,
                        required=True,
                        choices=("multi","individual","both"),
                        default="individual",
                        help="choose between having a multi-fa/gbk, invidual, or both for the output")

    parser.add_argument("--ret_type",
                        required=True,
                        choices=("fasta","genbank"),
                        help="return format of file")

    parser.add_argument("--sleepy",
                        type=int,
                        default=30,
                        help="Amount to delay a query to NCBI by")

    #parser.add_argument("--data",
    #                    type=lambda x: is_dir(parser,x,"results"),
    #                    default="results/data_accs.txt")

    """
    parser.add_argument("--multi_output",
                        type=argparse.FileType("w+"),
                        default="multi")
    """
    parser.add_argument("--galaxy_on",
                        action="store_true",
                        help="user to run galaxy like outputs")


    parser.add_argument("--data_name",
                        type=str,
                        default="data_accs.txt",
                        help="name of acc file")

    args = parser.parse_args()
    #print(args)
    # Write individual records
    #if not os.path.exists("results"):
        #os.mkdir("results")
    print(os.getcwd())
    if not os.path.exists("results"):
        os.mkdir("results")
    path = os.path.join("results",args.data_name)

    #with open(path,"w+") as f:
    #    f.writelines("accessions: "+str(args.input)+"\n")

    #if args.galaxy_on:
    #    os.chdir("results")

    if "__at__" in args.email:
        splits = args.email.split("__at__")
        email = splits[0]+"@"+splits[1]
    elif "@" in args.email:
        email = args.email
    elif args.email is None:
        raise Exception("EMAIL IS NECESSARY TO USE TOOL")

    #  Join together admin emails to append to hopefully catch NCBI's eye if abuse occurs
    admins = ["curtisross@tamu.edu","cory.maughmer@tamu.edu","anthonyc@tamu.edu"]
    sep = ";"
    admins.insert(0,email)
    emails = sep.join(admins)

    print("Logged in as: "+email)
    count = 0 # add a counter, so, it will do a two minute delay every 20th query, to attempt to not bother NCBI with load.
    path = os.path.join("results","output")
    for acc in args.input:
        count += 1
        if count % 20 == 0:
            sleep(120)
            pass
        else:
            pass
        c = CPTEfetch(emails, acc, args.db, args.ret_type)
        print(c)
        if args.galaxy_on:
            c.write_record(st=args.sleepy,name=path,galaxy=True)
        else:
            c.write_record(st=args.sleepy,name="data_",galaxy=False)

    # If more multi format is requested, perform below
    if args.ret_format == "multi" or args.ret_format == "both":
        if args.galaxy_on:
            #awk_files("DAT",output=f"outputMulti.{str(args.ret_type)}")
            #awk_files(str(args.ret_type),output=f"outputMulti.{str(args.ret_type)}")
            awk_files(str(args.ret_type),output=path,galaxy=True)
        else:
            awk_files(str(args.ret_type),output="outputMulti"+str(args.ret_type))
    print("---finish---")
    print(os.getcwd())

    ### Test for subdirs
    folders = []
    for r, d, f in os.walk(os.getcwd()):
        for folder in d:
            folders.append(os.path.join(r, folder))
    
    for fold in folders:
        print(fold)
    
    resultsfiles = [f for f in os.listdir("results") if os.path.isfile(os.path.join("results",f))]
    print("files in results: " +str(resultsfiles))
