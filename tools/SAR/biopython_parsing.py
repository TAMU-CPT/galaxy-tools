#!/usr/bin/env python
# Biopython parsing module. Uses in conjunction with the sar_finder script, and potential future scripts down the line.

from Bio import SeqIO

class FASTA_parser:
    """ Parses multi fasta file, and zips together header with sequence """

    def __init__(self, fa):
        self.fa = fa
    
    def multifasta_dict(self):
        """ parses the input multi fasta, and puts results into dictionary """

        return SeqIO.to_dict(SeqIO.parse(self.fa,"fasta"))


if __name__ == "__main__":
    fa_file = "test-data/mu-proteins.fa"
    d = FASTA_parser(fa_file).multifasta_dict()
    print(d)
    for k, v in d.items():
        print(v.description)
    
