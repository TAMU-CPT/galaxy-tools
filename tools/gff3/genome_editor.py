#!/usr/bin/env python
import logging
import copy
import argparse
import tsv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from CPT_GFFParser import gffParse, gffWrite
from gff3 import feature_lambda, feature_test_contains

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def mutate(gff3, fasta, changes, customSeqs, new_id):
    # Change Language
    # - we can only accept ONE genome as an input. (TODO: support multiple?)
    # - we can only build ONE genome as an output. (TODO: support multiple?)
    # - must allow selection of various regions
    # '1,1000,+   40,100,-    custom_seq_1'
    try:
        custom_seqs = SeqIO.to_dict(SeqIO.parse(customSeqs, "fasta"))
    except:
        custom_seqs = {}
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    # Pull first and onl record
    rec = list(gffParse(gff3, base_dict=seq_dict))[0]
    # Create a "clean" record
    new_record = copy.deepcopy(rec)
    new_record.id = new_id
    new_record.seq = Seq("")
    new_record.features = []
    new_record.annotations = {}
    # Process changes.
    chain = []
    topFeats = {}
    covered = 0
    for feat in rec.features:
        if "ID" in feat.qualifiers.keys():
          topFeats[feat.qualifiers["ID"][0]] = feat.location.start
    for change in changes:
        if "," in change:
            (start, end, strand) = change.split(",")
            start = int(start) - 1
            end = int(end)

            # Make any complaints
            broken_feature_start = list(
                feature_lambda(
                    rec.features,
                    feature_test_contains,
                    {"index": start},
                    subfeatures=False,
                )
            )
            if len(broken_feature_start) > 0:
                pass
                # log.info("DANGER: Start index chosen (%s) is in the middle of a feature (%s %s). This feature will disappear from the output", start, broken_feature_start[0].id, broken_feature_start[0].location)
            broken_feature_end = list(
                feature_lambda(
                    rec.features,
                    feature_test_contains,
                    {"index": end},
                    subfeatures=False,
                )
            )
            if len(broken_feature_end) > 0:
                pass
                # log.info("DANGER: End index chosen (%s) is in the middle of a feature (%s %s). This feature will disappear from the output", end, broken_feature_end[0].id, broken_feature_end[0].location)

            # Ok, fetch features
            if strand == "+":
                tmp_req = rec[start:end]
            else:
                tmp_req = rec[start:end].reverse_complement(
                    id=True,
                    name=True,
                    description=True,
                    features=True,
                    annotations=True,
                    letter_annotations=True,
                    dbxrefs=True,
                )

            def update_location(feature, shiftS):
                feature.location = FeatureLocation(feature.location.start + shiftS, feature.location.end + shiftS, feature.strand)
                for i in feature.sub_features:
                  i = update_location(i, shiftS)
                return feature
                
            

            #for feature in tmp_req.features:
            
                  
                

            chain.append(
                [
                    rec.id,
                    start + 1,
                    end,
                    strand,
                    new_record.id,
                    len(new_record) + 1,
                    len(new_record) + (end - start),
                    "+",
                ]
            )

            covered += len(new_record.seq)
            print(covered)
            new_record.seq += tmp_req.seq
            # NB: THIS MUST USE BIOPYTHON 1.67. 1.68 Removes access to
            # subfeatures, which means you will only get top-level features.
            startInd = len(new_record.features)
            new_record.features += tmp_req.features
            
            for i in new_record.features[startInd:]:
                i.location = FeatureLocation(i.location.start + covered, i.location.end + covered, i.location.strand)
                if "ID" not in i.qualifiers.keys():
                  continue
                diffS = i.location.start - topFeats[i.qualifiers["ID"][0]]
                subFeats = i.sub_features
                for j in subFeats:
                  j = update_location(j, diffS)
        else:
            new_record.seq += custom_seqs[change].seq
    yield new_record, chain


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=argparse.FileType("r"), help="Sequence")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="Annotations")
    parser.add_argument("new_id", help="Append to ID", default="_v2")
    parser.add_argument(
        "--out_fasta",
        type=argparse.FileType("w"),
        help="Output fasta",
        default="out.fa",
    )
    parser.add_argument(
        "--out_gff3",
        type=argparse.FileType("w"),
        help="Output gff3",
        default="out.gff3",
    )
    parser.add_argument(
        "--out_simpleChain",
        type=argparse.FileType("w"),
        help="Output simple chain (i.e. not a real UCSC chain file)",
        default="out.chain",
    )
    parser.add_argument("--changes", nargs="+")
    parser.add_argument("--customSeqs", type=argparse.FileType("r"))
    args = parser.parse_args()

    for rec, chain in mutate(
        args.gff3, args.fasta, args.changes, args.customSeqs, args.new_id
    ):
        # TODO: Check that this appends and doesn't overwirte
        gffWrite([rec], args.out_gff3)
        SeqIO.write([rec], args.out_fasta, "fasta")
        tsv.dump(chain, args.out_simpleChain)
