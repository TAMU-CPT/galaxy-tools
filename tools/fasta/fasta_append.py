#!/usr/bin/env python
######## Takes 2 input FASTAs, and pairs headers to headers and sequences to sequences
######## FOR CONTROLS USE: test-data/appendA.fa and appendB.fa
######## Currently this will truncate an unevenly sized FASTA
######## An example would be the following:
'''
FASTA_a + FASTA_b -->  Header_1a+Header_1b | sequence_1a + sequence_1b ==> FASTA_1c
'''
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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='append two fasta files into one')

    parser.add_argument('fa_1', type=argparse.FileType("r"), help='FASTA file #1') # the first "input" argumen
    parser.add_argument('fa_2', type=argparse.FileType("r"), help='FASTA file #2') # the first "input" argumen
    parser.add_argument('--appended_fa', dest='appended_fa', type=argparse.FileType("w"), default='appended.fa', help='appended FASTA output') # appended file
    args = parser.parse_args()
    
    fa_1 = args.fa_1 # fasta file 1
    fa_2 = args.fa_2 # fasta file 2

    fa1 = split_fa(fa_1) # zipped data of fasta 1 ; (description, sequence)
    fa2 = split_fa(fa_2) # zipped data of fasta 2 ; (description, sequence)

    fastas = zip(fa1,fa2) # zip the zip

    with args.appended_fa as f: # open file from args.appended_fa
        for fasta_pairs in fastas: # parse through the zipped zip
            f.write('>'+str(fasta_pairs[0][0])+'_'+str(fasta_pairs[1][0])) # fasta file 1 header_ fasta file 2 header
            f.write('\n'+lineWrapper(str(fasta_pairs[0][1]+fasta_pairs[1][1]))+'\n') # fasta file 1 seq + fasta file 2 seq



