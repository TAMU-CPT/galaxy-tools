#!/usr/bin/env python3
import sys
import time
import argparse

from apollo import accessible_organisms
from apollo.util import GuessOrg, OrgOrGuess
from arrow.apollo import get_apollo_instance

from webapollo import WAAuth, OrgOrGuess, GuessOrg, AssertUser, retry
from cpt_gffParser import gffParse, gffWrite
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sample script to add an attribute to a feature via web services"
    )
    parser.add_argument("apollo", help="Complete Apollo URL")
    parser.add_argument("username", help="WA Username")
    parser.add_argument("password", help="WA Password")
    parser.add_argument("email", help="User Email")
    parser.add_argument("--source", help="URL where the input dataset can be found.")
    parser.add_argument(
        "--name",
        default="product",
        help="Qualifer to use for the name of the added feature, default is 'product'",
    )
    parser.add_argument(
        "--org_json",
        type=argparse.FileType("r"),
        help="Apollo JSON output, source for common name",
    )
    parser.add_argument("--org_raw", help="Common Name")
    parser.add_argument("--org_id", help="Organism ID")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 file")
    
    args = parser.parse_args()

    #wa = WebApolloInstance(args.apollo, args.username, args.password)

    wa = get_apollo_instance()
    # User must have an account
    #gx_user = AssertUser(wa.users.loadUsers(email=args.email))

    # Get organism
    if args.org_json:
        orgs = [x.get("commonName", None) for x in json.load(args.org_json)]
        orgs = [x for x in orgs if x is not None]
        org_cn = orgs
    elif args.org_raw:
        org = args.org_raw.strip()
        if len(org) > 0:
            org_cn = [org]
        else:
            raise Exception("Organism Common Name not provided")
    elif args.org_id:
        org_cn = [wa.organisms.findOrganismById(args.org_id).get("commonName", None)]
    else:
        raise Exception("Organism Common Name not provided")

    annoteClient = wa.annotations

    if isinstance(org_cn, list):
        org_cn = org_cn[0]

    # TODO: Check user perms on org.
    org = wa.organisms.findOrganismByCn(org_cn)

    
#    (self, organism, gff3, 
    source=args.source
    batch_size=1
    test=False,
    use_name=False,
    disable_cds_recalculation=False,
    timing=False,
        """
        Load a full GFF3 into annotation track
        :type organism: str
        :param organism: Organism Common Name
        :type gff3: str
        :param gff3: GFF3 to load
        :type source: str
        :param source: URL where the input dataset can be found.
        :type batch_size: int
        :param batch_size: Size of batches before writing
        :type test: bool
        :param test: Run in dry run mode
        :type use_name: bool
        :param use_name: Use the given name instead of generating one.
        :type disable_cds_recalculation: bool
        :param disable_cds_recalculation: Disable CDS recalculation and instead use the one provided
        :type timing: bool
        :param timing: Output loading performance metrics
        :rtype: str
        :return: Loading report
        """
    organisms = wa.organisms.get_organisms()
    org_ids = []
    for org in organisms:
        if org_cn == org['commonName'] or org_cn == str(org['id']):
            org_ids.append(org['id'])

        if len(org_ids) == 0:
            raise Exception("Organism name or id not found [" + org_cn + "]")

        if len(org_ids) > 1:
            raise Exception("More than one organism found for [" + org_cn + "].  Use an organism ID instead: " + str(org_ids))

    total_features_written = 0
        
    all_processed = {'top-level': [], 'transcripts': []}
    loading_status = {}
    for rec in gffParse(args.gff3):
        annoteClient.set_sequence(organism, rec.id)
        try:
            log.info("Processing %s with features: %s" % (rec.id, rec.features))
            processed = annoteClient._process_gff_entry(rec, source=source,
                                                disable_cds_recalculation=disable_cds_recalculation,
                                                use_name=use_name
                                                )
            all_processed['top-level'].extend(processed['top-level'])
            all_processed['transcripts'].extend(processed['transcripts'])
            total_features_written += 1
            written_top = annoteClient._check_write(batch_size, test, all_processed['top-level'], FeatureType.FEATURE, timing)
            written_transcripts = annoteClient._check_write(batch_size, test, all_processed['transcripts'], FeatureType.TRANSCRIPT, timing)

            if len(written_top):
                all_processed['top-level'] = []
                loading_status = {**loading_status, **written_top}
            if len(written_transcripts):
                all_processed['transcripts'] = []
                loading_status = {**loading_status, **written_transcripts}

        except Exception as e:
            msg = str(e)
            if '\n' in msg:
                msg = msg[0:msg.index('\n')]
            log.error("Failed to load features from %s" % rec.id)

        # Write the rest of things to write (ignore batch_size)
    written_top = annoteClient._check_write(0, test, all_processed['top-level'], FeatureType.FEATURE, timing)
    written_transcripts = annoteClient._check_write(0, test, all_processed['transcripts'], FeatureType.TRANSCRIPT, timing)

    if len(written_top):
        all_processed['top-level'] = []
        loading_status = {**loading_status, **written_top}
    if len(written_transcripts):
        all_processed['transcripts'] = []
        loading_status = {**loading_status, **written_transcripts}

    log.info("Finished loading")
        

#    return loading_status
