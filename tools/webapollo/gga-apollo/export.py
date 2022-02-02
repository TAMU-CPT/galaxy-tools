#!/usr/bin/env python
from __future__ import print_function

import argparse
import json
import time

from apollo import accessible_organisms
from apollo.util import CnOrGuess, GuessCn, GuessOrg

from arrow.apollo import get_apollo_instance

from webapollo import UserObj, handle_credentials

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to export data from Apollo via web services')
    CnOrGuess(parser)
    parser.add_argument('--gff', type=argparse.FileType('w'))
    parser.add_argument('--seq', type=argparse.FileType('w'))
    parser.add_argument('--fasta_pep', type=argparse.FileType('w'))
    parser.add_argument('--fasta_cds', type=argparse.FileType('w'))
    parser.add_argument('--fasta_cdna', type=argparse.FileType('w'))
    parser.add_argument('--vcf', type=argparse.FileType('w'))
    parser.add_argument('--json', type=argparse.FileType('w'))
    parser.add_argument('email', help='User Email')
    args = parser.parse_args()

    wa = get_apollo_instance()

    # User must have an apollo account, if not, create it
    gx_user = UserObj(**wa.users._assert_or_create_user(args.email))
    handle_credentials(gx_user)

    org_cn = []
    if args.org_json:
        org_cn = [x.get("commonName", None) for x in json.load(args.org_json)]
        org_cn = [x for x in orgs if x is not None]
    elif args.org_raw:
        args.org_raw = str(args.org_raw)
        args.org_raw = args.org_raw.split(",")
        for i in range(len(args.org_raw)):
            org_cn.append(args.org_raw[i].strip())
    elif args.org_id:
        orgList = str(args.org_id)
        orgList = orgList.split(",")
        for x in orgList:
            args.org_id = x.strip()
            res = GuessOrg(args, wa)
            if res:
               org_cn.append(res[0])

    if len(org_cn) == 0:
        raise Exception("Organism Common Name not provided")


    org_data = []
    for org_ind in org_cn:
          if args.org_raw:
            args.org_raw = org_ind.strip()
          elif args.org_id:
            args.org_id = org_ind.strip()
          orgs, seqs = GuessCn(args, wa)    
          org = wa.organisms.show_organism(org_ind)

          if args.gff:
              uuid_gff = wa.io.write_downloadable(org['commonName'], 'GFF3', export_gff3_fasta=True, sequences=seqs)
              if 'error' in uuid_gff or 'uuid' not in uuid_gff:
                  raise Exception("Apollo failed to prepare the GFF3 file for download: %s" % uuid_gff)
              args.gff.write(wa.io.download(uuid_gff['uuid'], output_format="text"))
              time.sleep(1)

          if args.seq:
              uuid_fa = wa.io.write_downloadable(org['commonName'], 'FASTA', sequences=seqs, seq_type='genomic')
              if 'error' in uuid_fa or 'uuid' not in uuid_fa:
                  raise Exception("Apollo failed to prepare the FASTA file for download: %s" % uuid_fa)
              args.seq.write(wa.io.download(uuid_fa['uuid'], output_format="text"))
              time.sleep(1)

          if args.vcf:
              uuid_vcf = wa.io.write_downloadable(org['commonName'], 'VCF', sequences=seqs)
              if 'error' in uuid_vcf or 'uuid' not in uuid_vcf:
                  raise Exception("Apollo failed to prepare the VCF file for download: %s" % uuid_vcf)
              args.vcf.write(wa.io.download(uuid_vcf['uuid'], output_format="text"))
              time.sleep(1)

          if args.fasta_cdna:
              uuid_fa = wa.io.write_downloadable(org['commonName'], 'FASTA', sequences=seqs, seq_type='cdna')
              if 'error' in uuid_fa or 'uuid' not in uuid_fa:
                  raise Exception("Apollo failed to prepare the cdna FASTA file for download: %s" % uuid_fa)
              args.fasta_cdna.write(wa.io.download(uuid_fa['uuid'], output_format="text"))
              time.sleep(1)

          if args.fasta_cds:
              uuid_fa = wa.io.write_downloadable(org['commonName'], 'FASTA', sequences=seqs, seq_type='cds')
              if 'error' in uuid_fa or 'uuid' not in uuid_fa:
                  raise Exception("Apollo failed to prepare the cds FASTA file for download: %s" % uuid_fa)
              args.fasta_cds.write(wa.io.download(uuid_fa['uuid'], output_format="text"))
              time.sleep(1)

          if args.fasta_pep:
              uuid_fa = wa.io.write_downloadable(org['commonName'], 'FASTA', sequences=seqs, seq_type='peptide')
              if 'error' in uuid_fa or 'uuid' not in uuid_fa:
                  raise Exception("Apollo failed to prepare the file for download: %s" % uuid_fa)
              args.fasta_pep.write(wa.io.download(uuid_fa['uuid'], output_format="text"))

          org_data.append(org)

    args.json.write(json.dumps(org_data, indent=2))
