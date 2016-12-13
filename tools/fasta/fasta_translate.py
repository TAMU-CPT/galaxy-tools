#!/usr/bin/env python
import sys
import logging
import argparse
from Bio import SeqIO
from Bio.Data import CodonTable
logging.basicConfig(level=logging.INFO)


def translate(fasta_file, target='protein', table=11, strip_stops=False):
    records = list(SeqIO.parse(fasta_file, "fasta"))

    # If we strip stops, then OK as is, otherwise add the extra char that the
    # [0:] would've missed
    strip_stops_diff = 0 if strip_stops else -1

    for record in records:
        if target == 'protein':
            mod = len(record.seq) % 3
            if mod != 0:
                record.seq = record.seq[0:-mod]

            try:
                tmpseq = record.seq.translate(table=table, cds=True)
            except CodonTable.TranslationError:
                tmpseq = record.seq.translate(table=table, cds=False)

            if '*' in tmpseq:
                tmpseq = tmpseq[0:str(tmpseq).index('*') - strip_stops_diff]

            record.seq = tmpseq
            if len(record.seq) > 0:
                SeqIO.write(record, sys.stdout, "fasta")
        else:
            record.seq = record.seq.transcribe()
            SeqIO.write(record, sys.stdout, "fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Translate fasta file')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), help='Fasta file')
    parser.add_argument('--target', choices=['protein', 'rna'])
    parser.add_argument('--table', type=int, default=11,
                        help='Translation table to use', choices=range(1, 23))
    parser.add_argument('--strip_stops', action='store_true', help='Remove stop characters')

    args = parser.parse_args()
    translate(**vars(args))
