'''
PREMISE
### Functions/Classes that are used in both generate-putative-osp.py and generate-putative-isp.py
###### Main premise here is to make the above scripts a little more DRY, as well as easily readable for execution.
###### Documentation will ATTEMPT to be thourough here
'''

import re
from Bio import SeqIO
from Bio import Seq

# Not written in OOP for a LITTLE bit of trying to keep the complication down in case adjustments are needed by someone else.
# Much of the manipulation is string based; so it should be straightforward as well as moderately quick
################## GLobal Variables
Lys = 'K'

def check_back_end_snorkels(seq, tmsize):
    """
        Searches through the backend of a potential TMD snorkel. This is the 2nd part of a TMD snorkel lysine match.
        --> seq : should be the sequence fed from the "search_region" portion of the sequence
        --> tmsize : size of the potential TMD being investigated
    """
    found = []
    if (seq[tmsize-4] == Lys and re.search(('[FIWLVMYCATGS]'),seq[tmsize-5])):
        found = 'match'
        return found
    elif (seq[tmsize-3] == Lys and re.search(('[FIWLVMYCATGS]'),seq[tmsize-4])):
        found = 'match'
        return found
    elif (seq[tmsize-2] == Lys and re.search(('[FIWLVMYCATGS]'),seq[tmsize-3])):
        found = 'match'
        return found
    elif (seq[tmsize-1] == Lys and re.search(('[FIWLVMYCATGS]'),seq[tmsize-2])):
        found = 'match'
        return found
    else:
        found = 'NOTmatch'
        return found



def find_tmd(pair,minimum=10,maximum=30,TMDmin=10,TMDmax=20):
    """ 
        Function that takes a string input AA sequence and searches for TMDs:
        ---> seq : Input of tuple with description and AA sequence (str)
        ---> minimum : How close from the initial start codon a TMD can be within
        ---> maximum : How far from the initial start codon a TMD can be within
        Function that searches for lysine snorkels and then for a spanning hydrophobic region that indicates a potential TMD
        ---> TMDmin : The minimum size that a transmembrane can be (default = 10)
        ---> TMDmax : The maximum size tha ta transmembrane can be (default = 20)
    """
    # hydrophobicAAs = ['F', 'I', 'W', 'L', 'V', 'M', 'Y', 'C', 'A', 'T', 'G', 'S']
    tmd = []
    s = str(pair[1]) # sequence being analyzed
    #print(s) # for trouble shooting
    search_region = s[minimum-1:maximum+1]
    #print(search_region) # for trouble shooting
    for tmsize in range(TMDmin, TMDmax+1, 1):
        #print('==============='+str(tmsize)+'================') # print for troubleshooting
        pattern = "['FIWLVMYCATGS']{"+str(tmsize)+"}" # searches for these hydrophobic residues tmsize total times
        if re.search(('[K]'), search_region[0:7]):
            try:
                backend = check_back_end_snorkels(search_region,tmsize)
                if backend == 'match':
                    tmd.append(pair)
                else:
                    continue
            except (IndexError,TypeError):
                continue
        elif re.search((pattern), search_region):
            try:
                tmd.append(pair)
            except (IndexError,TypeError):
                continue
        else:
            continue

        return tmd

def find_lipobox(pair,minimum=10,maximum=30,regex=1):
    """
        #### Function that takes an input tuple, and will return pairs of sequences to their description that have a lipoobox
        ####### minimum - min distance from start codon to first AA of lipobox
        ####### maximum - max distance from start codon to first AA of lipobox
        ####### regex - option 1 (default) => more strict regular expression ; option 2 => looser selection, imported from LipoRy
        
    """
    if regex == 1:
        pattern = '[ILMFTV][^REXD][GAS]C' # regex for Lipobox from findSpanin.pl
    elif regex == 2:
        pattern = '[ACGSILMFTV][^REKD][GASNL]C' # regex for Lipobox from LipoRy

    candidates = []
    s = str(pair[1])
    #print(s) # trouble shooting
    search_region = s[minimum:maximum+1]
    print(search_region) # trouble shooting
    #for each_pair in pair:
    print(s)
    if re.search((pattern), search_region): # lipobox must be WITHIN the range...
        # searches the sequence with the input RegEx AND omits if start sequence begins with I
        candidates.append(pair)
        print('passed') # trouble shooting
        return candidates
    else:
        print('didnotpass') # trouble shooting
        pass

def tuple_fasta(fasta_file):
    """
        #### INPUT: Fasta File
        #### OUTPUT: zipped (zip) : pairwise relationship of description to sequence
        #### 
    """
    fasta = SeqIO.parse(fasta_file, 'fasta')
    descriptions = []
    sequences = []
    for r in fasta: # iterates and stores each description and sequence
        description = (r.description) 
        sequence = str(r.seq)
        if sequence[0] != 'I': # the translation table currently has I as a potential start codon ==> this will remove all ORFs that start with I
            descriptions.append(description)
            sequences.append(sequence)
        else:
            continue
    
    return zip(descriptions, sequences)
    

####################################################################################################################################
if __name__ == "__main__":

#### TMD TEST
    test_desc = ['one','two','three','four','five']
    test_seq = ['XXXXXXXXXXXXXXXFMCFMCFMCFMCFMCXXXXXXXXXXXXXXXXXXXXXXXXXX','XXXXXXXXAAKKKKKKKKKKKKKKKXXXXXXXXXXXXX','XXXXXXX','XXXXXXXXXXXKXXXXXXXXXX','XXXXXXXXXXAKXXXXXXXXXXAKXXXXXXXX']
    #for l in 
    #combo = zip(test_desc,test_seq)
    pairs = zip(test_desc, test_seq)
    tmd = []
    for each_pair in pairs:
        print(each_pair)
        try:
            tmd += find_tmd(pair=each_pair)
        except (IndexError,TypeError):
            continue
        #try:s = each_pair[1]
        #tmd += find_tmd(seq=s, tmsize=15)
    print('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    print(tmd)
    print('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    
    

    #### tuple-fasta TEST
        #fasta_file = 'out_isp.fa'
        #ret = tuple_fasta(fasta_file)
        #print('=============')
        #for i in ret:
            #print(i[1])

    #### LipoBox TEST
    test_desc = ['one','two','three','four','five','six','seven']
    test_seq = ['XXXXXXXXXTGGCXXXXXXXXXXXXXXXX','XXXXXXXXAAKKKKKKKKKKKKKKKXXXXXXXXXXXXX','XXXXXXX','AGGCXXXXXXXXXXXXXXXXXXXXTT','XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXTGGC','XXXXXXXXXXXXXXXXXXXXXXXXXXTGGC','MSTLRELRLRRALKEQSMRYLLSIKKTLPRWKGALIGLFLICVATISGCASESKLPEPPMVSVDSSLMVEPNLTTEMLNVFSQ*']
    pairs = zip(test_desc, test_seq)
    lipo = []
    for each_pair in pairs:
        print(each_pair)
        #try:
        try:
            lipo += find_lipobox(pair=each_pair, regex=2)#, minimum=8)
        except TypeError: # catches if something doesnt have the min/max requirements (something is too small)
            continue
        #except:
            #continue
    print('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    print(lipo)
    print('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
