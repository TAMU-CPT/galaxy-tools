'''
PREMISE
### Functions/Classes that are used in both generate-putative-osp.py and generate-putative-isp.py
###### Main premise here is to make the above scripts a little more DRY, as well as easily readable for execution.
###### Documentation will ATTEMPT to be thourough here
'''

import re
from Bio import SeqIO
from Bio import Seq

def find_tmd(pair,start=10,stop=30,tmsize=16):
    """ 
        Function that takes a string input AA sequence and searches for TMDs:
        ---> seq : Input of tuple with description and AA sequence (str)
        ---> Start : How close from the initial start codon a TMD can be within
        ---> End : How far from the initial start codon a TMD can be within
    """
    Lys = 'K'*3
    # hydrophobicAAs = ['F', 'I', 'W', 'L', 'V', 'M', 'Y', 'C', 'A', 'T', 'G', 'S']
    pattern = "['FIWLVMYCATGS']{"+str(tmsize)+"}" # searches for these hydrophobic residues tmsize total times
    tmd = []
    s = str(pair[1])
    if re.search((pattern), s[start:stop+1]):
        tmd.append(pair)
    else: # re.search((pattern), seq[start:stop+tmsize+1]):
        sequence = pair[1][start:stop+tmsize+1]
        if sequence[0:2] == Lys:
            tmd.append(pair)
        elif sequence[tmsize-6:tmsize-3] == Lys:
            tmd.append(pair)
        else:
            pass
    return tmd

def find_lipobox(pair,start=10,stop=30):
    """
        #### Function that takes an input tuple, and will return pairs of sequences to their description that have a lipoobox
        ####### start - min distance from start codon to first AA of lipobox
        ####### stop - max distance from start codon to first AA of lipobox
        
    """
    regggie = '[ACGSILMFTV][^REKD][GASNL]C' # regex for Lipobox from LipoRy
    pattern = '[ILMFTV][^REXD][GAS]C' # regex for Lipobox from findSpanin.pl

    candidates = []
    for each_pair in pairs:
        s = str(each_pair[1])
        if (re.search((pattern), s[start:stop]) and s[0] != 'I'):
            # searches the sequence with the input RegEx AND omits if start sequence begins with I
            candidates.append(each_pair)
        else:
            continue

    return candidates

def tuple_fasta(fasta_file):
    """
        #### INPUT: Fasta File
        #### OUTPUT: zipped (zip) of each pairwise description to sequence
    """
    fasta = SeqIO.parse(fasta_file, 'fasta')
    descriptions = []
    sequences = []
    for r in fasta: # iterates and stores each description and sequence
        description = (r.description) 
        sequence = str(r.seq)
        if sequence[0] != 'I':
            descriptions.append(description)
            sequences.append(sequence)
        else:
            continue
    
    return zip(descriptions, sequences)
    

####################################################################################################################################
if __name__ == "__main__":

    test_desc = ['one','two','three']
    test_seq = ['XXXXXXXXXXXXXXXFMCFMCFMCFMCFMCXXXXXXXXXXXXXXXXXXXXXXXXXX','XXXXXXXXKKKKKKKKKKKKKKKKXXXXXXXXXXXXX','XXXXXXX']
    #combo = zip(test_desc,test_seq)
    pairs = zip(test_desc, test_seq)
    tmd = []
    for each_pair in pairs:
        tmd += find_tmd(pair=each_pair, tmsize=15 )
        #try:s = each_pair[1]
        #tmd += find_tmd(seq=s, tmsize=15)
    #print(tmd)
    
    #print('++++++++++++')

    fasta_file = 'out_isp.fa'
    ret = tuple_fasta(fasta_file)
    #print('=============')
    #for i in ret:
        #print(i[1])