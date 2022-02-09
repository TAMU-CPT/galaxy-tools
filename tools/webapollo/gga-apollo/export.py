#!/usr/bin/env python
from __future__ import print_function

import argparse
import json
import time

from apollo import accessible_organisms
from apollo.util import CnOrGuess, GuessCn, GuessOrg

from arrow.apollo import get_apollo_instance

from webapollo import UserObj, handle_credentials


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to export data from Apollo via web services"
    )
    CnOrGuess(parser)
    parser.add_argument('--gff', type=argparse.FileType('w'))
    parser.add_argument('--gff_with_fasta', action='store_true')
    parser.add_argument("--seq", type=argparse.FileType("w"))
    parser.add_argument("--fasta_pep", type=argparse.FileType("w"))
    parser.add_argument("--fasta_cds", type=argparse.FileType("w"))
    parser.add_argument("--fasta_cdna", type=argparse.FileType("w"))
    parser.add_argument("--vcf", type=argparse.FileType("w"))
    parser.add_argument("--json", type=argparse.FileType("w"))
    parser.add_argument('--die', action='store_true')
    parser.add_argument("email", help="User Email")
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

    if len(org_cn) > 5:
        # Don't want to be abused by a user that has lots of organisms
        # and that decides to export/delete them all in one go.
        # Make them spread it out.
        raise Exception("Please reduce organism count to 5 or less. Thank you!")

    def error(message):
        if args.die:
            raise Exception(message)
        else:
            print(message)

    org_data = []
    for org_ind in org_cn:
        if args.org_raw:
            args.org_raw = org_ind.strip()
        elif args.org_id:
            args.org_id = org_ind.strip()
        orgs, seqs = GuessCn(args, wa) #  this seqs is usually plugged into the write_downloadable sequence arg
        orgs_check = accessible_organisms(gx_user, [org_ind], 'READ')
        if not orgs_check:
            raise Exception("You do not have read permission on organism %s" % org_cn)
        org = wa.organisms.show_organism(org_ind)


        # Fetch all of the refseqs
        realSeqs = wa.organisms.get_sequences(org["id"])

        for sequence in realSeqs['sequences']:
            print("Downloading", sequence)

            try:
              if args.gff:
                uuid_gff = wa.io.write_downloadable(org['commonName'], 'GFF3', export_gff3_fasta=args.gff_with_fasta, sequences=[sequence['name']])
                if 'error' in uuid_gff or 'uuid' not in uuid_gff:
                    error("Apollo failed to prepare the GFF3 file for download: %s" % uuid_gff)
                args.gff.write(wa.io.download(uuid_gff['uuid'], output_format="text"))
                time.sleep(1)
            except Exception as e:
              error(e)

            try:
              if args.vcf:
                uuid_vcf = wa.io.write_downloadable(org['commonName'], 'VCF', sequences=[sequence['name']])
                if 'error' in uuid_vcf or 'uuid' not in uuid_vcf:
                    error("Apollo failed to prepare the VCF file for download: %s" % uuid_vcf)
                args.vcf.write(wa.io.download(uuid_vcf['uuid'], output_format="text"))
                time.sleep(1)
            except Exception as e:
              error(e)

            try:
              if args.fasta_cdna:
                uuid_fa = wa.io.write_downloadable(org['commonName'], 'FASTA', sequences=[sequence['name']], seq_type='cdna')
                if 'error' in uuid_fa or 'uuid' not in uuid_fa:
                    error("Apollo failed to prepare the cdna FASTA file for download: %s" % uuid_fa)
                args.fasta_cdna.write(wa.io.download(uuid_fa['uuid'], output_format="text"))
                time.sleep(1)
            except Exception as e:
              error(e)

            try:
              if args.fasta_cds:
                uuid_fa = wa.io.write_downloadable(org['commonName'], 'FASTA', sequences=[sequence['name']], seq_type='cds')
                if 'error' in uuid_fa or 'uuid' not in uuid_fa:
                    error("Apollo failed to prepare the cds FASTA file for download: %s" % uuid_fa)
                args.fasta_cds.write(wa.io.download(uuid_fa['uuid'], output_format="text"))
                time.sleep(1)
            except Exception as e:
              error(e)

            try:
              if args.fasta_pep:
                uuid_fa = wa.io.write_downloadable(org['commonName'], 'FASTA', sequences=[sequence['name']], seq_type='peptide')
                if 'error' in uuid_fa or 'uuid' not in uuid_fa:
                    error("Apollo failed to prepare the file for download: %s" % uuid_fa)
                args.fasta_pep.write(wa.io.download(uuid_fa['uuid'], output_format="text"))
                time.sleep(1)
            except Exception as e:
              error(e)

            org_data.append(org)

    args.json.write(json.dumps(org_data, indent=2))
