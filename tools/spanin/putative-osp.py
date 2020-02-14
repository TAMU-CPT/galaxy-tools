import argparse
from cpt import OrfFinder
from Bio import SeqIO
from Bio import Seq
from BCBio import GFF
from lipory import find_lipoprotein
import os
import sys

### Requirement Inputs
#### INPUT : Genomic FASTA
#### OUTPUT : Putative OSP candidates in GFF3 and FASTA format.
### Notes:
####### As of 2.13.2020 - RegEx pattern: ^.{%s,%s}[ACGSILMFTV][^REKD][GASNL]C for Lipobox
if __name__ == '__main__':

    # Common parameters for both ISP / OSP portion of scripts

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
    
    # osp parameters
    parser.add_argument('--osp_min_len', dest='osp_min_len', default=30, help='Minimum ORF length, measured in codons', type=int)
    parser.add_argument('--osp_on', dest='out_osp_nuc', type=argparse.FileType('w'), default='out_osp.fna', help='Output nucleotide sequences, FASTA')
    parser.add_argument('--osp_op', dest='out_osp_prot', type=argparse.FileType('w'), default='out_osp.fa', help='Output protein sequences, FASTA')
    parser.add_argument('--osp_ob', dest='out_osp_bed', type=argparse.FileType('w'), default='out_osp.bed', help='Output BED file')
    parser.add_argument('--osp_og', dest='out_osp_gff3', type=argparse.FileType('w'), default='out_osp.gff3', help='Output GFF3 file')
    parser.add_argument('--osp_min_dist', dest='osp_min_dist', default=10, help='Minimal distance to first AA of TMD, measured in AA', type=int)
    parser.add_argument('--osp_max_dist', dest='osp_max_dist', default=30, help='Maximum distance to first AA of TMD, measured in AA', type=int)
    parser.add_argument('--putative_osp', dest='putative_osp_gff3', type=argparse.FileType('w'), default='putative_osp.gff3', help='Output of putative GFF3 file')

    parser.add_argument('-v', action='version', version='0.3.0') # Is this manually updated?
    args = parser.parse_args()
    the_args = vars(parser.parse_args())
    
    ### osp output, naive ORF finding:
    osps = OrfFinder(args.table, args.ftype, args.ends, args.osp_min_len, args.strand)
    osps.locate(args.fasta_file, args.out_osp_nuc, args.out_osp_prot, args.out_osp_bed, args.out_osp_gff3)

    print('++++name++++')
    print(args.fasta_file.name)
    print('^^^^name^^^^')
    for k, v in the_args.items():
        print(k+' : '+str(v))
    
    slimwithlipo = find_lipoprotein(args.out_osp_gff3.name, 
                                    args.fasta_file.name, 
                                    lipobox_mindist=args.osp_min_dist, 
                                    lipobox_maxdist=args.osp_max_dist, 
                                    spanin=True) # Chopping overall noise with LipoRy

    for record in slimwithlipo:
        record[0].annotations = {}
        #print(record)
        GFF.write(record, sys.stdout)
