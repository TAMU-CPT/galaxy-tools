import argparse
import sys
from time import sleep
from Bio import Entrez
from Bio import SeqIO
import urllib


# Not very DRY, but, I think this will get the job done and be easy to read/alter for anyone in the future.
def do_multi(user_input,db,sleep_amt,ret_type):
    """ Formats a multi-fasta/gbk """
    acc_list = user_input
    try:
        if ret_type == "genbank":
            net_handle = Entrez.efetch(db=db,id=acc_list,rettype="gb",retmode="text")
        else:
            net_handle = Entrez.efetch(db=db,id=acc_list,rettype="fasta",retmode="text")
    except urllib.error.HTTPError: # <-- This is the response error if you get timed out, if it occurs, sleep and send the request again
        sleep(sleep_amt)
        if ret_type == "genbank":
            net_handle = Entrez.efetch(db=db,id=acc_list,rettype="gb",retmode="text")
        else:
            net_handle = Entrez.efetch(db=db,id=acc_list,rettype="fasta",retmode="text")
    records = SeqIO.parse(net_handle,ret_type)
    with open(f"multi{ret_type}.{ret_type}", "w") as file:
        for record in records:
            SeqIO.write(record, file, ret_type)

    return print("...File(s) fetched and downloaded...")

def do_individuals(user_input,db,sleep_amt,ret_type):
    """ returns individual fasta/gbk files """
    acc_list = user_input
    for acc in acc_list:
        try:
            if ret_type == "genbank":
                net_handle = Entrez.efetch(db=db,id=acc,rettype="gbwithparts", retmode="text")
            else:
                net_handle = Entrez.efetch(db=db,id=acc,rettype=ret_type, retmode="text")
        except urllib.error.HTTPError: # <-- This is the response error if you get timed out, if it occurs, sleep and send the request again
            sleep(sleep_amt)
            if ret_type == "genbank":
                net_handle = Entrez.efetch(db=db,id=acc,rettype="gbwithparts", retmode="text")
            else:
                net_handle = Entrez.efetch(db=db,id=acc,rettype=ret_type, retmode="text")
        record = net_handle.read()
        print(record)
        with open(f"{acc}.{ret_type}","w") as file:
            file.write(record)
        net_handle.close()
        file.close()

    return print("...File(s) fetched and downloaded...")




if __name__ == "__main__":

    ##### Arguments
    parser = argparse.ArgumentParser(description="CPT's very own modified Efetch")

    parser.add_argument("email",
                        type=str,
                        help="Entrez Required Email") # current place holder until I determine how best to use the current user's email from Galaxy

    parser.add_argument("--input",
                        type=str,
                        nargs="*",
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
    Entrez.email = args.email
    if args.ret_format == "multi":
        do_multi(user_input=args.input,db=args.db,ret_type=args.ret_type,sleep_amt=args.sleep)
    elif args.ret_format == "individual":
        do_individuals(user_input=args.input,db=args.db,ret_type=args.ret_type,sleep_amt=args.sleep)
    elif args.ret_format == "both":
        do_multi(user_input=args.input,db=args.db,ret_type=args.ret_type,sleep_amt=args.sleep)
        do_individuals(user_input=args.input,db=args.db,ret_type=args.ret_type,sleep_amt=args.sleep)