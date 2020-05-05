import time
from Bio import Entrez
from Bio import SeqIO
from urllib.error import HTTPError
import argparse
from helperFunctions import cat_files


Entrez.email = "curtisross@tamu.edu"

class CPTEfetch:
    """ Object that has built in functions to retrieve data from NCBI. Initially constructued to retreive GB and FA files from the nuccore and protein NCBI databases """

    def __init__(self,email,acc,db,ret_type):
        self.email = email
        self.acc = acc
        self.db = db
        self.ret_type = ret_type


    def __repr__(self):
        return "<accession: {}| database: {}>".format(self.acc,self.db)


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


    def write_record(self,st,galaxy=True):
        record = self.retrieve_data(sleep_time=st)
        if galaxy:
            with open(f"{str(self.acc)}.DAT","w") as file:
                file.write(record)
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
                        action="append",
                        #nargs="*",
                        help='accession input"')

    parser.add_argument("--db",
                        type=str,
                        choices=("protein", "nuccore"),
                        help="choose protein or nuccore database to do query")

    parser.add_argument("--ret_format",
                        type=str,
                        choices=("multi","individual","both"),
                        default="individual",
                        help="choose between having a multi-fa/gbk, invidual, or both for the output")

    parser.add_argument("--ret_type",
                        choices=("fasta","genbank"),
                        help="return format of file")

    parser.add_argument("--sleep",
                        type=int,
                        default=20,
                        help="Amount to delay a query to NCBI by")

    parser.add_argument("--output",
                        type=argparse.FileType("w+"),
                        default="output.dat")

    args = parser.parse_args()

    print(args.input)
    for i in args.input:
        print(i)
    print(args.email)
    print(args.db)
    print(args.ret_type)

    for acc in args.input:
        c = CPTEfetch(args.email, acc, args.db, args.ret_type)
        #c = CPTEfetch("curtisross@tamu.edu","NC_000866.4","nuccore","genbank")
        print(c)
        c.write_record(st=args.sleep,galaxy=False)

    if args.ret_format == "multi" or args.ret_format == "both":
        cat_files(args.ret_type)