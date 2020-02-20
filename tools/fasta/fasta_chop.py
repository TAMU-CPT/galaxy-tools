################################################################################
###### Goal : Takes a FASTA file and trims sequence to a given value ###########
################################################################################
################## Use chop.fa for a control to test. (fasta/test-data/chop.fa)

from Bio import SeqIO
from Bio import Seq
import argparse

def split_fa(fa): # recycled from spaninFuncs.py
    """
        Splits a fasta into a tuple; (description, sequence)
        ---> INPUT : fasta file
        ---> OUTPUT : tuple of fasta separation
    """
    fasta = SeqIO.parse(fa, 'fasta')
    descriptions = []
    sequences = []
    for r in fasta: # iterates and stores each description and sequence
        description = (r.description) 
        sequence = str(r.seq)
        descriptions.append(description)
        sequences.append(sequence)    
    return zip(descriptions, sequences)

def lineWrapper(text, charactersize=60): # recycled from spaninFuncs.py
    """
        wraps text for file inputs. handy for fasta files.
    """
    if len(text) <= charactersize:
        return text
    else:
        return text[:charactersize] + '\n' + lineWrapper(text[charactersize:], charactersize)

def chop_it(fa_tup, amt=30):
    """
        Takes input tuple zip (from split fa) and chops the sequences to a given size
    """
    chopped = []
    for each_pair in fa_tup:
        seq = each_pair[1][:amt]
        chop = (each_pair[0],seq)
        chopped.append(chop)
    return chopped

if __name__ == "__main__":
    ##### Arugments to be fed:
    parser = argparse.ArgumentParser(description='append two fasta files into one')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='fasta file') # the first "input" argument
    parser.add_argument('--chopped_fa', dest='chopped_fa', type=argparse.FileType("w"), default='chopped.fa', help='chopped fasta output') # appended file
    parser.add_argument('--chop_amt', dest='chop_amt', type=int, default=30, help='amount to chop a fasta file by')
    args = parser.parse_args()

    fasta = args.fasta

    tupes = split_fa(fasta)

    chopped = chop_it(tupes, amt=args.chop_amt)

    with args.chopped_fa as f: # open file from args.appended_fa
        for data in chopped: # parse through the zipped zip
            f.write('>'+str(data[0])) # fasta file 1 header_ fasta file 2 header
            f.write('\n'+str(data[1])+'\n') # fasta file 1 seq + fasta file 2 seq
