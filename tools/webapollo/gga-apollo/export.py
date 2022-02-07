#!/usr/bin/env python
from __future__ import print_function

import argparse
import json
import time

from apollo import accessible_organisms
from apollo.util import CnOrGuess, GuessCn, GuessOrg

from arrow.apollo import get_apollo_instance

from webapollo import UserObj, handle_credentials


def check_file_can_be_written(obj, seq_type, export_type):
    """
    Checks if a file can be written based on if the returned object from Apollo is a null object or not.

    If the object is null, it will exit the conditional block with an empty file.

    This is to allow the delete organism tool to still process. TODO: We should modify it to write out in the file *why* it is empty.

    If another error occurs, it will raise the exception and error out as before.

    :type obj: IOClient object(?) TODO: Confirm
    :param obj: WebApollo IOClient object.

    :type seq_type: str
    :param seq_type: Export selection. Choices: peptide, cds, cdna, genomic, or empty (ie '')

    :type export_type: str
    :param export_type: Export type. Choices: FASTA, GFF3, VCF
    """

    matches = [
        "null object",
        #"index out of range",
    ]

    if any(match in obj["error"] for match in matches):
        pass
    else:
        raise Exception(
            f"Apollo failed to prepare the {seq_type} {export_type} file for download: {obj}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to export data from Apollo via web services"
    )
    CnOrGuess(parser)
    parser.add_argument("--gff", type=argparse.FileType("w"))
    parser.add_argument("--seq", type=argparse.FileType("w"))
    parser.add_argument("--fasta_pep", type=argparse.FileType("w"))
    parser.add_argument("--fasta_cds", type=argparse.FileType("w"))
    parser.add_argument("--fasta_cdna", type=argparse.FileType("w"))
    parser.add_argument("--vcf", type=argparse.FileType("w"))
    parser.add_argument("--json", type=argparse.FileType("w"))
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
        # and that decides to delete them all in one go.
        # Make them spread it out.
        raise Exception("Please reduce organism count to 5 or less. Thank you!")

    org_data = []
    for org_ind in org_cn:
        if args.org_raw:
            args.org_raw = org_ind.strip()
        elif args.org_id:
            args.org_id = org_ind.strip()
        orgs, seqs = GuessCn(args, wa)
        org = wa.organisms.show_organism(org_ind)

        if args.gff:
            uuid_gff = wa.io.write_downloadable(
                org["commonName"], "GFF3", export_gff3_fasta=True, sequences=seqs
            )
            if "error" in uuid_gff or "uuid" not in uuid_gff:
                check_file_can_be_written(uuid_gff, "", "GFF3")
            else:
                args.gff.write(wa.io.download(uuid_gff["uuid"], output_format="text"))
                time.sleep(1)

        if args.seq:
            uuid_fa = wa.io.write_downloadable(
                org["commonName"], "FASTA", sequences=seqs, seq_type="genomic"
            )
            if "error" in uuid_fa or "uuid" not in uuid_fa:
                check_file_can_be_written(uuid_fa, "genomic", "FASTA")
            else:
                args.seq.write(wa.io.download(uuid_fa["uuid"], output_format="text"))
                time.sleep(1)

        if args.vcf:
            uuid_vcf = wa.io.write_downloadable(
                org["commonName"], "VCF", sequences=seqs
            )
            if "error" in uuid_vcf or "uuid" not in uuid_vcf:
                check_file_can_be_written(uuid_vcf, "", "VCF")
            else:
                args.vcf.write(wa.io.download(uuid_vcf["uuid"], output_format="text"))
                time.sleep(1)

        if args.fasta_cdna:
            uuid_fa = wa.io.write_downloadable(
                org["commonName"], "FASTA", sequences=seqs, seq_type="cdna"
            )
            if "error" in uuid_fa or "uuid" not in uuid_fa:
                check_file_can_be_written(uuid_fa, "cdna", "FASTA")
            else:
                args.fasta_cdna.write(
                    wa.io.download(uuid_fa["uuid"], output_format="text")
                )
                time.sleep(1)

        if args.fasta_cds:
            uuid_fa = wa.io.write_downloadable(
                org["commonName"], "FASTA", sequences=seqs, seq_type="cds"
            )
            if "error" in uuid_fa or "uuid" not in uuid_fa:
                check_file_can_be_written(uuid_fa, "cds", "FASTA")
            else:
                args.fasta_cds.write(
                    wa.io.download(uuid_fa["uuid"], output_format="text")
                )
                time.sleep(1)

        if args.fasta_pep:
            uuid_fa = wa.io.write_downloadable(
                org["commonName"], "FASTA", sequences=seqs, seq_type="peptide"
            )
            if "error" in uuid_fa or "uuid" not in uuid_fa:
                check_file_can_be_written(uuid_fa, "peptide", "FASTA")
            else:
                args.fasta_pep.write(
                    wa.io.download(uuid_fa["uuid"], output_format="text")
                )

        org_data.append(org)

    args.json.write(json.dumps(org_data, indent=2))
