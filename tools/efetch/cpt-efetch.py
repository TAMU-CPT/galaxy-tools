##### 

import argparse
from time import sleep
from Bio import Entrez
from Bio import SeqIO


def query(email,user_input,db,sleep_amt,multi=False,both=False):
    """ Use Entrez.efetch to query based on a user's request """

    Entrez.email = email
    print("Logged in as: "+str(Entrez.email))

    if multi:
        acc_list = user_input
        net_handle = Entrez.efetch(db=db,id=acc_list,rettype="fasta",retmode="txt")
    else:
        acc_list = user_input
        for acc in acc_list:
            try: # sleep this
                net_handle = Entrez.efetch(db=db,id=acc,rettype="fasta", retmode="text")
            except HTTPError:
                sleep(sleep_amt)
                net_handle = Entrez.efetch(db=db,id=acc,rettype="fasta",retmode="text")
    print(net_handle.read())
    sleep(sleep_amt)




if __name__ == "__main__":

    ##### Arguments
    parser = argparse.ArgumentParser(description="Trim the putative protein candidates and find potential i-spanin / o-spanin pairs")

    parser.add_argument("email",
                        type=str,
                        help="Entrez Required Email") # current place holder until I determine how best to use the current user's email from Galaxy

    parser.add_argument("--input",
                        type=str,
                        nargs="*",
                        help='accession input"')

    parser.add_argument("--multi",
                        action="store_true",
                        help="store as multifasta or multigbk")

    parser.add_argument("--both",
                        action="store_true",
                        help="do both single and multi fetching")

    parser.add_argument("--db",
                        type=str,
                        choices=("protein", "nuccore"),
                        help="choose protein or nuccore database to do query")
    
    parser.add_argument("--sleep",
                        type=int,
                        default=1,
                        help="Amount to delay a query to NCBI by")

    args = parser.parse_args()

    ##### NCBI Query
    query(email=args.email,user_input=args.input,db=args.db,sleep_amt=args.sleep)
