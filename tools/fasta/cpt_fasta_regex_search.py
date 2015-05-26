#!/usr/bin/env python
import re
import logging
logging.basicConfig(level=logging.INFO)
import argparse
from Bio import SeqIO
import StringIO

def search(fasta, regex):
    m = re.compile(regex)
    good_records = []
    for record in SeqIO.parse(fasta, "fasta"):
        if m.search(str(record.seq).upper()):
            good_records.append(record)
    return good_records

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter sequences based on regular expressions')
    parser.add_argument('fasta', type=file, help='Fasta Sequence')
    parser.add_argument('regex', help='Regular Expression')

    args = parser.parse_args()

    output = StringIO.StringIO()
    records = search(**vars(args))
    SeqIO.write(records, output, "fasta")
    print output.getvalue()
