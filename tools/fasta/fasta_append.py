
######## Takes X amount of input FASTAs, and pairs headers to headers and sequences to sequences
######## An example would be the following:
'''
FASTA_a + FASTA_b -->  Header_1a+Header_1b | sequence_1a + sequence_1b ==> FASTA_1c
'''
import re
from Bio import SeqIO
from Bio import Seq

def split_fa(fa):
    """
        Splits a fasta into a tuple; (description, sequence)
        ---> INPUT : fasta file
        ---> OUTPUT : tuple of fasta separation
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


