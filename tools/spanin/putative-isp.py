##### findSpanin.pl --> findSpanin.py
######### Incooperated from the findSpanin.pl script, but better and more snakey.

import argparse
from cpt import OrfFinder
from Bio import SeqIO
from Bio import Seq
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

    parser.add_argument('--table', dest='table', default=1, 
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

    parser.add_argument('-v', action='version', version='0.3.0') # Is this manually updated?
    args = parser.parse_args()
    the_args = vars(parser.parse_args())
    
    ### isp output, naive ORF finding:
    isps = OrfFinder(args.table, args.ftype, args.ends, args.isp_min_len, args.strand)
    isps.locate(args.fasta_file, args.out_isp_nuc, args.out_isp_prot, args.out_isp_bed, args.out_isp_gff3)

    print('++++name++++')
    print(args.fasta_file.name)
    print('^^^^name^^^^')
    for k, v in the_args.items():
        print(k+' : '+str(v))