#!/usr/bin/env python
import sys
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def search(fasta, regex):
    m = re.compile(regex)
    for record in SeqIO.parse(fasta, "fasta"):
        if m.search(str(record.seq).upper()):
            oseq = str(record.seq).upper()
            newseq = ""
            last_idx = 0
            for hit in m.finditer(oseq):
                newseq += oseq[last_idx:hit.start()].lower()
                newseq += oseq[hit.start():hit.end()].upper()
                last_idx = hit.end()
            newseq += oseq[last_idx:].lower()
            # Update sequence
            record.seq = Seq(newseq)
            yield [record]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter sequences based on regular expressions')
    parser.add_argument('fasta', type=file, help='Fasta Sequence')
    parser.add_argument('regex', help='Regular Expression (use upper case sequences)')

    args = parser.parse_args()

    for record in search(**vars(args)):
        SeqIO.write(record, sys.stdout, "fasta")
