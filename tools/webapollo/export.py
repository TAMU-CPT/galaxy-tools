#!/usr/bin/env python
import sys
import StringIO
import json
import argparse
from Bio import SeqIO
from BCBio import GFF
from webapollo import WAAuth, WebApolloInstance, CnOrGuess, GuessCn


def export(org_cn, seqs):
    org_data = wa.organisms.findOrganismByCn(org_cn)

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
        ).encode('utf-8'))
    else:
        data.write(wa.io.write(
            exportAllSequences=True,
            sequences=[],
            **kwargs
        ).encode('utf-8'))

    # Seek back to start
    data.seek(0)

    records = list(GFF.parse(data))
    if len(records) == 0:
        print "Could not find any sequences or annotations for this organism + reference sequence"
        sys.exit(2)
    else:
        for record in records:
            record.annotations = {}
            if args.gff:
                GFF.write([record], args.gff)
            record.description = ""
            if args.fasta:
                SeqIO.write([record], args.fasta, 'fasta')

    return org_data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    WAAuth(parser)
    CnOrGuess(parser)
    parser.add_argument('--gff', type=argparse.FileType('w'))
    parser.add_argument('--fasta', type=argparse.FileType('w'))
    parser.add_argument('--json', type=argparse.FileType('w'))

    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    org_cn_list, seqs = GuessCn(args)

    org_data = []
    for org_cn in org_cn_list:
        indiv_org_data = export(org_cn, seqs)
        org_data.append(indiv_org_data)
    args.json.write(json.dumps(org_data, indent=2))
