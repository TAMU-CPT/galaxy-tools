import time
import os
from Bio import Entrez
from Bio import SeqIO
from urllib.error import HTTPError
import argparse
from helperFunctions import awk_files


Entrez.email = "curtisross@tamu.edu"

class CPTEfetch:
    """ Object that has built in functions to retrieve data from NCBI. Initially constructued to retreive GB and FA files from the nuccore and protein NCBI databases """

    def __init__(self,email,acc,db,ret_type):
        self.email = email
        self.acc = acc
        self.db = db
        self.ret_type = ret_type


    def __repr__(self):
        return "<accession: {} | database: {} | return_type: {}>".format(self.acc,self.db,self.ret_type)


    def retrieve_data(self,sleep_time=5):
        if self.ret_type == "genbank":
            try:
                net_handle = Entrez.efetch(db=self.db,id=self.acc,rettype="gbwithparts",retmode="text")
            except HTTPError:
                time.sleep(sleep_time)
                net_handle = Entrez.efetch(db=self.db,id=self.acc,rettype="gbwithparts",retmode="text")
        else:
            try:
                net_handle = Entrez.efetch(db=self.db,id=self.acc,rettype="fasta",retmode="text")
            except HTTPError:
                time.sleep(sleep_time)
                net_handle = Entrez.efetch(db=self.db,id=self.acc,rettype="fasta",retmode="text")
        return net_handle.read()


    def write_record(self,name,st,galaxy=True):
        record = self.retrieve_data(sleep_time=st)
        if galaxy:
            with open(f"{name.name}_{str(self.acc)}.{str(self.ret_type)}","w") as file:
                file.write(record)
            #with open(f"{name}_{str(self.acc)}.DAT","w") as file:
                #file.write(record)
        else:
            with open(f"{str(self.acc)}.{str(self.ret_type)}","w") as file:
                file.write(record)


if __name__ == "__main__":
        ##### Arguments
    parser = argparse.ArgumentParser(description="CPT's very own modified Efetch")

    parser.add_argument("email",
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

    parser.add_argument("--sleep",
                        type=int,
                        default=20,
                        help="Amount to delay a query to NCBI by")

    parser.add_argument("--data",
                        type=argparse.FileType("w"),
                        default="output")

    """
    parser.add_argument("--multi_output",
                        type=argparse.FileType("w+"),
                        default="multi")
    """
    parser.add_argument("--galaxy_on",
                        action="store_true",
                        help="user to run galaxy like outputs")


    args = parser.parse_args()

    # Write individual records
    if args.galaxy_on:
        os.chdir("results")

    for acc in args.input:
        c = CPTEfetch(args.email, acc, args.db, args.ret_type)
        print(c)
        if args.galaxy_on:
            c.write_record(st=args.sleep,name=args.data,galaxy=True)
        else:
            c.write_record(st=args.sleep,name="results/output",galaxy=False)

    # If more multi format is requested, perform below
    if args.ret_format == "multi" or args.ret_format == "both":
        if args.galaxy_on:
            #awk_files("DAT",output=f"outputMulti.{str(args.ret_type)}")
            #awk_files(str(args.ret_type),output=f"outputMulti.{str(args.ret_type)}")
            awk_files(str(args.ret_type),output=args.data,galaxy=True)
        else:
            awk_files(str(args.ret_type),output=f"outputMulti.{str(args.ret_type)}")
