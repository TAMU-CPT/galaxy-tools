#!/usr/bin/env python
import StringIO
import sys
import json
import argparse
from Bio import SeqIO
from BCBio import GFF
from webapollo import WAAuth, WebApolloInstance

if __name__ == '__main__':
    json
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)

    parser.add_argument('commonName', help='Common Name')
    parser.add_argument('refSeqs', nargs='*', help='SeqName(s)')

    parser.add_argument('--gff', type=argparse.FileType('w'))
    parser.add_argument('--fasta', type=argparse.FileType('w'))

    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    data = StringIO.StringIO()

    kwargs = dict(
        exportType='GFF3',
        seqType='genomic',
        exportGff3Fasta=True,
        output="text",
        exportFormat="text",
        organism=args.commonName,
    )

    for refSeq in args.refSeqs:
        data.write(wa.io.write(
            exportAllSequences=False,
            sequences=refSeq,
            **kwargs
        ))
    else:
        data.write(wa.io.write(
            exportAllSequences=True,
            sequences=[],
            **kwargs
        ))

    # Seek back to start
    data.seek(0)

    for record in GFF.parse(data):
        record.annotations = {}
        GFF.write([record], args.gff)
        record.description = ""
        SeqIO.write([record], args.fasta, 'fasta')
