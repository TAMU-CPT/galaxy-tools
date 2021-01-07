#!/usr/bin/env python
import argparse
import json
import logging

from apollo import util, accessible_organisms
from apollo.client import Client
from apollo.util import features_to_feature_schema, retry, GuessOrg, OrgOrGuess
from apollo.annotations import FeatureType
from enum import Enum

from arrow.apollo import get_apollo_instance

from cpt_gffParser import gffParse

from webapollo import UserObj, handle_credentials
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

#class FeatureType(Enum):
#    FEATURE = 1
#    TRANSCRIPT = 2

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('email', help='User Email')
    parser.add_argument('--source', help='URL where the input dataset can be found.')
    parser.add_argument('--use_name', action='store_true', help='Use the given name instead of generating one.')
    parser.add_argument('--disable_cds_recalculation', action='store_true', help='Disable CDS recalculation and instead use the one provided.')
    parser.add_argument('--overrideID', type=str, help='Select new field to display name in annotation track')
    OrgOrGuess(parser)

    parser.add_argument('gff3', type=argparse.FileType('r'), help='GFF3 file')
    args = parser.parse_args()

    wa = get_apollo_instance()
    # User must have an account
    gx_user = UserObj(**wa.users._assert_or_create_user(args.email))
    handle_credentials(gx_user)

    # Get organism
    org_cn = GuessOrg(args, wa)
    if isinstance(org_cn, list):
        org_cn = org_cn[0]

    all_orgs = wa.organisms.get_organisms()
    if 'error' in all_orgs:
        all_orgs = []
    all_orgs = [org['commonName'] for org in all_orgs]
    if org_cn not in all_orgs:
        raise Exception("Could not find organism %s" % org_cn)

    orgs = accessible_organisms(gx_user, [org_cn], 'WRITE')
    if not orgs:
        raise Exception("You do not have write permission on this organism")

    #load_result = wa.annotations.load_legacy_gff3(org_cn, args.gff3) # use_name=args.use_name, disable_cds_recalculation=args.disable_cds_recalculation)
    #print(json.dumps(load_result, indent=2))

    annoteClient = wa.annotations
    all_processed = {'top-level': [], 'transcripts': []}
    total_features_written = 0
    batch_size = 1
    test = False
    timing=False
    loading_status = {}
    
    recList = [] 
    for rec in gffParse(args.gff3):
        filteredFeats = []  
        featList = rec.features
        featList.sort(key=lambda x: x.location.start)
        for x in featList:
            for y in ["pseudogene", "pseudogenic_region", "processed_pseudogene", 'transcript', 'tRNA', 'snRNA', 'snoRNA', 'ncRNA', 'rRNA', 'mRNA', 'miRNA', 'guide_RNA', 'RNase_P_RNA', 'telomerase_RNA', 'SRP_RNA', 'lnc_RNA', 'RNase_MRP_RNA', 'scRNA', 'piRNA', 'tmRNA', 'enzymatic_RNA', "repeat_region", "terminator", "shine_dalgarno_sequence", "transposable_element", "gene", "CDS"]:
                if str(x.type) == y:
                    filteredFeats.append(x)
                    if args.overrideID != "ID" and args.overrideID in filteredFeats[-1].qualifiers.keys():
                        filteredFeats[-1].qualifiers["ID"] = filteredFeats[-1].qualifiers[args.overrideID]
                        filteredFeats[-1].id = filteredFeats[-1].qualifiers["ID"][0]
                    break
        rec.features = filteredFeats
        recList.append(rec)

    #args.gff3.seek(0)
    for rec in recList:
        annoteClient.set_sequence(org_cn, rec.id)
        try:
            log.info("Processing %s with features: %s" % (rec.id, rec.features))
            processed = annoteClient._process_gff_entry(rec, source=args.source,
                                                disable_cds_recalculation=args.disable_cds_recalculation,
                                                use_name=args.use_name
                                                )
            print(processed)
            all_processed['top-level'].extend(processed['top-level'])
            all_processed['transcripts'].extend(processed['transcripts'])
            total_features_written += 1
            written_top = annoteClient._check_write(batch_size, test, all_processed['top-level'])#, FeatureType.FEATURE, timing)
            written_transcripts = annoteClient._check_write(batch_size, test, all_processed['transcripts'], FeatureType.TRANSCRIPT)#, timing)

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
    written_top = annoteClient._check_write(0, test, all_processed['top-level'])#, FeatureType.FEATURE, timing)
    written_transcripts = annoteClient._check_write(0, test, all_processed['transcripts'])#, FeatureType.TRANSCRIPT, timing)

    if len(written_top):
        all_processed['top-level'] = []
        loading_status = {**loading_status, **written_top}
    if len(written_transcripts):
        all_processed['transcripts'] = []
        loading_status = {**loading_status, **written_transcripts}

    log.info("Finished loading")
    print(loading_status)
