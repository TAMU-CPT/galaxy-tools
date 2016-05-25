#!/usr/bin/env python
import StringIO
import sys
import json
import argparse
from Bio import SeqIO
from BCBio import GFF
from webapollo import WAAuth, WebApolloInstance, CnOrGuess, GuessCn

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)
    CnOrGuess(parser)
    parser.add_argument('--gff', type=argparse.FileType('w'))
    parser.add_argument('--fasta', type=argparse.FileType('w'))
    parser.add_argument('--json', type=argparse.FileType('w'))

    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    org_cn, seqs = GuessCn(args)
    org_data = wa.organisms.findOrganismByCn(org_cn)
    args.json.write(json.dumps([org_data], indent=2))

    data = StringIO.StringIO()

    kwargs = dict(
        exportType='GFF3',
        seqType='genomic',
        exportGff3Fasta=True,
        output="text",
        exportFormat="text",
        organism=org_cn,
    )

    if len(seqs) > 0:
        data.write(wa.io.write(
            exportAllSequences=False,
            sequences=seqs,
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
        if args.gff:
            GFF.write([record], args.gff)
        record.description = ""
        if args.fasta:
            SeqIO.write([record], args.fasta, 'fasta')
