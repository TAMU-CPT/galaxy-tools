# Takes an input data structure and file, to determine if items form the datastructure exist or not.

import re
import argparse
import Bio
import argparse
from BCBio import GFF
import explodeJSON as ej

#### Text Files
class FileParser:
    """
       Parses a file, using different methods based on the _TYPE_ of file that it is
    """
    def __init__(self, file):
        self.file = file  # term that wants to be queried


#### GFF3
    def read_a_GFF3(self):
        '''
        Class for parsing a GFF3 and finding proximity terms
        '''
        gff = open(self.file)
        gff_recs = []
        for rec in GFF.parse(gff):
            print(rec)
            gff_recs.append(rec)
        gff.close()


    #### Genbank
    def read_GBK(self):
        pass


    #### FASTA
    def read_FASTA(self):
        pass

class Terms:
    def __init__(self, termType):
        self.termType = termType
    
    def formQuery(self,options=[],filename=None,custom=None,combo=False):
        """ based on the termType, the formation of a query list is made """
        if combo == False:
            if self.termType == "dbase":
                db_path = "data/lysis-family-expanded.json"
                db = ej.explodeJSON(db_path)
                db = db.readJSON()
                # Since this is currently a static dbase, I'm going to hardcode the key options that users have. In the future, there might be some shenanigans done to enhance/improve this choice
                if options != []:
                    query = options
                    terms = []
                    for q in query:
                        print(q)
                        terms.extend(db[q]) 
                else:
                    terms = []
                    for vals in db.values():
                        terms.extend(vals)
                
                print(terms)
            elif self.termType == "custom_text":
                # I don't know the best way to test this before wrapping. However, I'll take a stab at what I think it will be.
                terms = []
                terms.extend(custom)
            elif self.termType == "file":
                # read in a new line separate file.
                terms = open(filename).read().splitlines()
                print(terms)
        else:
             
            pass

        return terms
    
if __name__ == "__main__":
    # PARAMS
    parser = argparse.ArgumentParser(description="file location")
    parser.add_argument("--termType",default="dbase")
    parser.add_argument("--options",default=[],nargs="*")
    parser.add_argument("--custom",default=None)
    parser.add_argument("--filename",default=None)
    args = parser.parse_args()
    
    # Function Calls
    t = Terms(termType=args.termType) # triggers response
    t.formQuery(options=args.options,custom=args.custom,filename=args.filename)