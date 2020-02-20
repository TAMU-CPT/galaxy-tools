#!/usr/bin/env python
import sys
try:
    import StringIO as io
except ImportError:
    import io

import json
import argparse
from Bio import SeqIO
from BCBio import GFF
##import gffutils
from webapollo import WAAuth, WebApolloInstance, CnOrGuess, GuessCn


def export(org_cn, seqs):
    org_data = wa.organisms.findOrganismByCn(org_cn)

    data = io.StringIO()

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
    #print(type(data))
    #print(dir(data))
    #print(data.getvalue())
    #exit()

    if args.gff:
      mode = 0
    else if args.fasta:
      mode = -1
    else:
      return org_data

    line = data.readline()
    while line:
      if line[0:2] == '..':
        if args.fasta:
          mode = 1
        else:
          return org_data
        line = data.readline()
        args.fasta.write('>' + line + '\n')
        line = data.readLine()
      elif (line [0:3] == '###'):
        line = data.readLine() # continue
      elif mode == 0:          
        args.gff.write(line + '\n')
        line = data.readLine()
      elif mode == 1:
        args.fasta.write(line + '\n')
        line = data.readLine()
      elif mode == -1:
        line = data.readLine()
      else:
        print("Unaccounted for line: " + line)
        line = data.readLine()

    #records = list(GFF.parse(data))
##    db = gffutils.create_db(data, dbfn='temp.db', force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
##    db2 = gffutils.FeatureDB('temp.db', keep_order=True)
    if False == True:
        print("Could not find any sequences or annotations for this organism + reference sequence")
        sys.exit(2)
##    else:
        #for record in records:
        #    record.annotations = {}
        #    record.features = sorted(record.features, key=lambda x: x.location.start)
##            if args.gff:
##              for y in db2.directives:
##                args.gff.write('##' + y +'\n')
##              for x in db2.all_features():
##                args.gff.write(x + '\n')
                #args.gff.writeGFF.write([record], args.gff)
            #record.description = ""
##            if args.fasta:
                
##                args.fasta.write([record], args.fasta, 'fasta')

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

    org_cn_list, seqs = GuessCn(args, wa)

    org_data = []
    for org_cn in org_cn_list:
        indiv_org_data = export(org_cn, seqs)
        org_data.append(indiv_org_data)
    args.json.write(json.dumps(org_data, indent=2))
