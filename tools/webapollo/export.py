#!/usr/bin/env python
import StringIO
import sys
import json
import argparse
from Bio import SeqIO
from BCBio import GFF
from webapollo import WebApolloInstance
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

if __name__ == '__main__':
    json
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Username')
    parser.add_argument('password', help='WA Password')

    parser.add_argument('commonName', nargs='+', help='Sequence Unique Names')

    parser.add_argument('--gff', type=argparse.FileType('a'), default='out.gff3')
    parser.add_argument('--fasta', type=argparse.FileType('a'), default='out.fa')

    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    for cn in args.commonName:
        log.info("Fetching %s" % cn)
        data = StringIO.StringIO(wa.io.write(
            exportType='GFF3',
            seqType='genomic',
            exportAllSequences=False,
            exportGff3Fasta=True,
            output="text",
            exportFormat="text",
            # TODO: CPT specific convention!!!!!!!!
            organism=[cn],
            sequences=[cn],
        ).encode('utf-8'))
        data.seek(0)

        for record in GFF.parse(data):
            record.annotations = {}
            GFF.write([record], args.gff)
            record.description = ""
            SeqIO.write([record], args.fasta, 'fasta')
