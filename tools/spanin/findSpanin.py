##### findSpanin.pl --> findSpanin.py
######### Incooperated from the findSpanin.pl script, but better and more snakey.

import argparse
import regex as re
from cpt import OrfFinder
from Bio import SeqIO
from Bio import Seq
import os

### Requirement Inputs
#### INPUT : Genomic FASTA
#### PARAMETERS :
######## strand : +, -, both
######## start codons: ATG, GTG and TTG
######## isp_min: minimal length of the ORF, measured in AAs
######## isp_nterm_mindist: minimal distance to first AA of TMD, measured in AA
######## isp_nterm_maxdist: maximum distance to first AA of TMD, measured in AA
######## osp_min: minimal length of the ORF, measured in AAs
######## osp_signal_mindist: minimal distance to first AA of Lipobox, measured in AA
######## osp_signal_maxdist: maximum distance to first AA of Lipobox, measured in AA
######## Use LipoRy (?)
######## max_isp_osp_distance: maximum distance between the END of the isp, 
########                       and the beginning of the osp, measured in AA

###############################################################################
############:::::::::::::        PART I           :::::::::::::################
###############::::: Parse the FASTA for potential ORFs :::::::################
###############################################################################

###############################################################################
############:::::::::::::        PART II           :::::::::::::###############
###############:::::::        Find best candidate       :::::::################
###############################################################################

###############################################################################
############:::::::::::::        PART III           :::::::::::::##############
###############::::: Output FASTA files, and best candidate :::::::############
###############################################################################

#if __name__ == '__main__':
    #pass
###############################################################################

if __name__ == '__main__':

    print(os.getcwd())

    parser = argparse.ArgumentParser(description='Get open reading frames')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), help='Fasta file')

    parser.add_argument('-f', '--format', dest='seq_format', default='fasta',
                        help='Sequence format (e.g. fasta, fastq, sff)')
    parser.add_argument('--table', dest='table', default=1, help='NCBI Translation table',
                        type=int)
    parser.add_argument('-t', '--ftype', dest='ftype', choices=('CDS', 'ORF'), default='ORF', help='Find ORF or CDSs')
    parser.add_argument('-e', '--ends', dest='ends', choices=('open', 'closed'), default='closed',
                        help='Open or closed. Closed ensures start/stop codons are present')
    parser.add_argument('-m', '--mode', dest='mode', choices=('all', 'top', 'one'), default='all',
                        help='Output all ORFs/CDSs from sequence, all ORFs/CDSs '
                        'with max length, or first with maximum length')
    parser.add_argument('--min_len', dest='min_len', default=10, help='Minimum ORF/CDS length', type=int)
    parser.add_argument('--strand', dest='strand', choices=('both', 'forward', 'reverse'), default='both', help='select strand')

    parser.add_argument('--on', dest='out_nuc', type=argparse.FileType('w'), default='out.fna', help='Output nucleotide sequences')
    parser.add_argument('--op', dest='out_prot', type=argparse.FileType('w'), default='out.fa', help='Output protein sequences')
    parser.add_argument('--ob', dest='out_bed', type=argparse.FileType('w'), default='out.bed', help='Output BED file')
    parser.add_argument('--og', dest='out_gff3', type=argparse.FileType('w'), default='out.gff3', help='Output GFF3 file')
    parser.add_argument('-v', action='version', version='0.3.0')
    args = parser.parse_args()

    of = OrfFinder(args.table, args.ftype, args.ends, args.min_len, args.strand)
    of.locate(args.fasta_file, args.out_nuc, args.out_prot, args.out_bed, args.out_gff3)










