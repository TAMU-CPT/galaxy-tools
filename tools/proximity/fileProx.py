# Takes an input data structure and file, to determine if items form the datastructure exist or not.

import re
import argparse
import Bio
from BCBio import GFF

#### Text Files
class FileParser:
    def __init__(self, file):
        self.file = file  # term that wants to be queried


#### GFF3
class read_a_GFF3(FileParser):
    def find_prox(self):
        gff = open(self.file)
        gff_recs = []
        for rec in GFF.parse(gff):
            print(rec)
            gff_recs.append(rec)
        gff.close()


#### Genbank
class read_a_GBK(FileParser):
    pass


#### FASTA
class read_a_FASTA(FileParser):
    pass


if __name__ == "__main__":
    pass
    ### What kind of input file are we going to read?

    ###### If it's a genbank, output neighboring results

    ###### If it's a gff3, output neighboring results

    ###### If it's a fast file, output neighboring results
