#!/usr/bin/env python
import json
import argparse
from webapollo import WebApolloInstance

if __name__ == '__main__':
    json
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Username')
    parser.add_argument('password', help='WA Password')

    parser.add_argument('commonName', nargs='+', help='Sequence Unique Names')

    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    # {"type":"GFF3", "seqType":"genomic", "exportAllSequences":"false",
     # "exportGff3Fasta":"true", "output":"file", "format":"gzip",
     # "sequences":["Merlin"]}

    print wa.io.write(
        exportType='GFF3',
        seqType='genomic',
        exportAllSequences=False,
        exportGff3Fasta=True,
        output="text",
        exportFormat="text",
        sequences=['Merlin']
    )



    # print json.dumps(wa.organisms.getSequencesForOrganism(
        # args.commonName
    # ), indent=4)
    # wa.annotations.setSequence("", "Bobbert")
    # print wa.annotations.getSequence(args.commonName, type='FASTA')
