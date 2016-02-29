#!/usr/bin/env python
import sys
import re
import argparse
import logging
logging.basicConfig()
log = logging.getLogger()
from Bio.Seq import Seq, reverse_complement, translate
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Data import CodonTable
from cpt import OrfFinder

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get open reading frames')
    parser.add_argument('fasta_file', type=file, help='Fasta file')

    parser.add_argument('-f', '--format', dest='seq_format',
                    default='fasta', help='Sequence format (e.g. fasta, fastq, sff)')
    parser.add_argument('--table', dest='table',
                    default=1, help='NCBI Translation table', type=int)
    parser.add_argument('-t', '--ftype', dest='ftype',
                    choices=('CDS', 'ORF'), default='ORF',
                    help='Find ORF or CDSs')
    parser.add_argument('-e', '--ends', dest='ends',
                    choices=('open', 'closed'), default='closed',
                    help='Open or closed. Closed ensures start/stop codons are present')
    parser.add_argument('-m', '--mode', dest='mode',
                    choices=('all', 'top', 'one'), default='all',
                    help='Output all ORFs/CDSs from sequence, all ORFs/CDSs '
                    'with max length, or first with maximum length')
    parser.add_argument('--min_len', dest='min_len',
                    default=10, help='Minimum ORF/CDS length', type=int)

    parser.add_argument('--on', dest='out_nuc',type=argparse.FileType('w'),
                    default='out.fna', help='Output nucleotide sequences')
    parser.add_argument('--op', dest='out_prot',type=argparse.FileType('w'),
                    default='out.fa', help='Output protein sequences',)
    parser.add_argument('--ob', dest='out_bed',type=argparse.FileType('w'),
                    default='out.bed', help='Output BED file')
    parser.add_argument('--og', dest='out_gff3', type=argparse.FileType('w'),
                    default='out.gff3', help='Output GFF3 file')
    parser.add_argument('-v', action='version', version='0.3.0')
    args = parser.parse_args()

    of = OrfFinder(args.table, argd.ftype, args.ends, args.min_len)
    of.locate(args.fasta_file, args.out_nuc, args.out_prot, args.out_bed, args.out_gff3)
