#!/usr/bin/env python
import StringIO
import sys
import json
import argparse
from Bio import SeqIO
from BCBio import GFF
from webapollo import WebApolloInstance

if __name__ == '__main__':
    json
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Username')
    parser.add_argument('password', help='WA Password')

    parser.add_argument('commonName', nargs='+', help='Sequence Unique Names')

    parser.add_argument('--gff', type=argparse.FileType('w'))
    parser.add_argument('--fasta', type=argparse.FileType('w'))

    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    data = StringIO.StringIO(wa.io.write(
        exportType='GFF3',
        seqType='genomic',
        exportAllSequences=False,
        exportGff3Fasta=True,
        output="text",
        exportFormat="text",
        sequences=args.commonName
    ))
    data.seek(0)

    for record in GFF.parse(data):
        record.annotations = {}
        GFF.write([record], args.gff)
        record.description = ""
        SeqIO.write([record], args.fasta, 'fasta')
        sys.exit()
