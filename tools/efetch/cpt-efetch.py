import time
from Bio import Entrez
from Bio import SeqIO
from urllib.error import HTTPError


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
        record = net_handle.read()
        #print(record)
        with open(f"{str(self.acc)}.{str(self.ret_type)}","w") as file:
            file.write(record)
        net_handle.close()  # probably redundant

if __name__ == "__main__":
    c = CPTEfetch("curtisross@tamu.edu","NC_001416.1","nuccore","fasta")
    print(c)
    c.retrieve_data()