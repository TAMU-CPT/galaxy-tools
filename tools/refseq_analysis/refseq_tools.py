import time
import sys
import os
from progress.bar import Bar

# External
from Bio import Entrez
from Bio import SeqIO

"""
The goal of this script is to insert a genome accession and determine 
if it has a RefSeq Record. Combinations of esearch/elink/esummary will be used from the Entrez suite.
esearch -db nuccore -query 'KC821634.1 FR687252.1' | elink -related -name nuccore_nuccore_gbrs | esummary | xtract -pattern DocumentSummary -element Extra,AssemblyAcc
"""

class CPTrefseq:
    
    def __init__(self, email, acc):

        self.email = email
        self.acc = acc
        Entrez.email = self.email
    
    def check_if_ref(self):
        """
        Pipes together esearch, elink, and esummary retrieve confirmation of existing RefSeq records
        """
        #  Check if current accession is already a RefSeq record (XX_)
        if type(self.acc) is int:
            gids = [self.acc]
        elif self.acc[2] == "_": # make sure we don't already have 
            return "Is_RefSeq", "Not_Fetched"
        #  Check ESearch and return results
        else:
            time.sleep(1.05) #  important to not overload NCBI (currently)
            record = self.do_esearch(self.acc)
            gids = record['IdList']

        #  ELink results from ESearch
        for gid in gids:
            record = self.do_elink(gid)
            try:
                linksetdb = record[0]['LinkSetDb'][0]
            except IndexError:
                return "No_Successful_Link", "Not_Fetched"
            if linksetdb['LinkName'] == 'nuccore_nuccore_gbrs':
                ids = linksetdb['Link']
                #  ESummary, verify RefSeq existence
                if not ids: # already empty, don't check
                    pass
                else:
                    for id in ids: # potential problem area
                        record = self.do_esummary(id['Id'])
                        ref_acc = record[0]['Caption']
                        title = record[0]['Title']
                        return ref_acc, title
            else:
                return "No_RefSeq_Record", "Not_Fetched"


    @staticmethod
    def do_esearch(acc):
        vals = Entrez.read(Entrez.esearch(db="nuccore",term=acc))
        
        return vals

    @staticmethod
    def do_elink(gid):
        data = Entrez.read(Entrez.elink(dbfrom="nuccore",db="nuccore",id=gid))

        return data
    
    @staticmethod
    def do_esummary(id):
        data = Entrez.read(Entrez.esummary(db="nuccore",id=id))
        
        return data


