##### findSpanin.pl --> findSpanin.py
######### Incooperated from the findSpanin.pl script, but better and more snakey.

import argparse
from cpt import OrfFinder
from Bio import SeqIO
from Bio import Seq
import re
from spaninFuncs import find_tmd, tuple_fasta
import os

#if __name__ == '__main__':
    #pass
###############################################################################

if __name__ == '__main__':

    # Common parameters for both ISP / OSP portion of script

    parser = argparse.ArgumentParser(description='Get putative protein candidates for spanins')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), 
                        help='Fasta file') # the "input" argument

    parser.add_argument('-f', '--format', dest='seq_format', default='fasta',
                        help='Sequence format (e.g. fasta, fastq, sff)') # optional formats for input, currently just going to do ntFASTA

    parser.add_argument('--strand', dest='strand', choices=('both', 'forward', 'reverse'), default='both', 
                        help='select strand') # Selection of +, -, or both strands

    parser.add_argument('--table', dest='table', default=11, 
                        help='NCBI Translation table',type=int) # Uses "default" NCBI codon table. This should always (afaik) be what we want...

    parser.add_argument('-t', '--ftype', dest='ftype', choices=('CDS', 'ORF'), default='ORF', 
                        help='Find ORF or CDSs') # "functional type(?)" --> Finds ORF or CDS, for this we want just the ORF

    parser.add_argument('-e', '--ends', dest='ends', choices=('open', 'closed'), default='closed',
                        help='Open or closed. Closed ensures start/stop codons are present') # includes the start and stop codon

    parser.add_argument('-m', '--mode', dest='mode', choices=('all', 'top', 'one'), default='all', # I think we want this to JUST be all...nearly always
                        help='Output all ORFs/CDSs from sequence, all ORFs/CDSs with max length, or first with maximum length')

    # isp parameters
    parser.add_argument('--isp_min_len', dest='isp_min_len', default=60, help='Minimum ORF length, measured in codons', type=int)
    parser.add_argument('--isp_on', dest='out_isp_nuc', type=argparse.FileType('w'), default='out_isp.fna', help='Output nucleotide sequences, FASTA')
    parser.add_argument('--isp_op', dest='out_isp_prot', type=argparse.FileType('w'), default='out_isp.fa', help='Output protein sequences, FASTA')
    parser.add_argument('--isp_ob', dest='out_isp_bed', type=argparse.FileType('w'), default='out_isp.bed', help='Output BED file')
    parser.add_argument('--isp_og', dest='out_isp_gff3', type=argparse.FileType('w'), default='out_isp.gff3', help='Output GFF3 file')
    parser.add_argument('--osp_min_dist', dest='isp_min_dist', default=10, help='Minimal distance to first AA of TMD, measured in AA', type=int)
    parser.add_argument('--osp_max_dist', dest='isp_max_dist', default=30, help='Maximum distance to first AA of TMD, measured in AA', type=int)
    parser.add_argument('--putative_isp', dest='putative_isp_fa', type=argparse.FileType('w'), default='putative_isp.fa', help='Output of putative FASTA file')
    
    parser.add_argument('-v', action='version', version='0.3.0') # Is this manually updated?
    args = parser.parse_args()
    the_args = vars(parser.parse_args())
    
    ### isp output, naive ORF finding:
    isps = OrfFinder(args.table, args.ftype, args.ends, args.isp_min_len, args.strand)
    isps.locate(args.fasta_file, args.out_isp_nuc, args.out_isp_prot, args.out_isp_bed, args.out_isp_gff3)

    ### read in robust ORF from OrfFinder/locate output
    #fasta = SeqIO.parse(args.out_isp_prot.name, 'fasta')
    '''
    descriptions = []
    sequences = []
    for r in fasta: # iterates and stores each description and sequence
        description = (r.description) 
        sequence = (r.seq)
        descriptions.append(description)
        sequences.append(sequence)
    '''
    pairs = tuple_fasta(fasta_file=args.out_isp_prot.name)
    #pairs = zip(descriptions,sequences) # combine the descriptions with their corresponding sequence from the two list

    candidates = [] # empty candidates list to be passed through the user input criteria
    for each_pair in pairs: # grab transmembrane domains based off of the spaninFuncts module; which follows loosely off transmembrane.py
        try:
            candidates += find_tmd(pair=each_pair, minimum=args.isp_min_dist, maximum=args.isp_max_dist)
        except TypeError:
            continue
    
    candidate_dict = { k:v for k,v in candidates}
    with open(args.putative_isp_fa.name, 'w') as f:
        for desc,s in candidate_dict.items(): # description / sequence
            f.write('> '+str(desc))
            f.write('\n'+str(s)+'\n')


    '''https://docs.python.org/3.4/library/subprocess.html'''
    '''https://github.tamu.edu/CPT/Galaxy-Tools/blob/f0bf4a4b8e5124d4f3082d21b738dfaa8e1a3cf6/tools/phage/transmembrane.py'''