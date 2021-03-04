#!/usr/bin/env python
import sys
import argparse
import logging
import uuid
from cpt_gffParser import gffParse, gffWrite
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from gff3 import feature_lambda, feature_test_type, get_id

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def main(fasta, gff3, feature_filter=None, nodesc=False):

    if feature_filter == "nice_cds":
        from gff2gb import gff3_to_genbank as cpt_Gff2Gbk

        for rec in cpt_Gff2Gbk(gff3, fasta, 11):
            seenList = {}
            if rec.seq[0] == "?":
                sys.stderr.write("Error: No Fasta ID matches GFF ID '" + rec.id + "'\n")
                exit(1)
            for feat in sorted(rec.features, key=lambda x: x.location.start):
                if feat.type != "CDS":
                    continue

                ind = 0
                if (
                    str(
                        feat.qualifiers.get("locus_tag", get_id(feat)).replace(" ", "-")
                    )
                    in seenList.keys()
                ):
                    seenList[
                        str(
                            feat.qualifiers.get("locus_tag", get_id(feat)).replace(
                                " ", "-"
                            )
                        )
                    ] += 1
                    ind = seenList[
                        str(
                            feat.qualifiers.get("locus_tag", get_id(feat)).replace(
                                " ", "-"
                            )
                        )
                    ]
                else:
                    seenList[
                        str(
                            feat.qualifiers.get("locus_tag", get_id(feat)).replace(
                                " ", "-"
                            )
                        )
                    ] = 1
                append = ""
                if ind != 0:
                    append = "_" + str(ind)

                if nodesc:
                    description = ""
                else:
                    feat.qualifiers["ID"] = [feat._ID]
                    product = feat.qualifiers.get("product", "")
                    description = "{1} [Location={0.location};ID={0.qualifiers[ID][0]}]".format(
                        feat, product
                    )
                yield [
                    SeqRecord(
                        feat.extract(rec).seq,
                        id=str(
                            feat.qualifiers.get("locus_tag", get_id(feat)).replace(
                                " ", "-"
                            )
                        )
                        + append,
                        description=description,
                    )
                ]

    elif feature_filter == "unique_cds":
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        seen_ids = {}

        for rec in gffParse(gff3, base_dict=seq_dict):
            noMatch = True
            if "Alias" in rec.features[0].qualifiers.keys():
                lColumn = rec.features[0].qualifiers["Alias"][0]
            else:
                lColumn = ""
            for x in seq_dict:
                if x == rec.id or x == lColumn:
                    noMatch = False
            if noMatch:
                sys.stderr.write("Error: No Fasta ID matches GFF ID '" + rec.id + "'\n")
                exit(1)
            newfeats = []
            for feat in sorted(
                feature_lambda(
                    rec.features, feature_test_type, {"type": "CDS"}, subfeatures=False
                ),
                key=lambda f: f.location.start,
            ):
                nid = rec.id + "____" + feat.id
                if nid in seen_ids:
                    nid = nid + "__" + uuid.uuid4().hex
                feat.qualifiers["ID"] = [nid]
                newfeats.append(feat)
                seen_ids[nid] = True

                if nodesc:
                    description = ""
                else:
                    if feat.strand == -1:
                      important_data = {"Location": FeatureLocation(feat.location.start + 1, feat.location.end - feat.shift, feat.strand)}
                    else:
                      important_data = {"Location": FeatureLocation(feat.location.start + 1 + feat.shift, feat.location.end, feat.strand)}
                    if "Name" in feat.qualifiers:
                        important_data["Name"] = feat.qualifiers.get("Name", [""])[0]

                    description = "[{}]".format(
                        ";".join(
                            [
                                "{key}={value}".format(key=k, value=v)
                                for (k, v) in important_data.items()
                            ]
                        )
                    )
                #if feat.id == "CPT_Privateer_006.p01":
                #print(feat)
                #exit()
                
                if isinstance(feat.location, CompoundLocation):
                  finSeq = ""
                  if feat.strand == -1:
                    for x in feat.location.parts:
                      finSeq += str((rec.seq[feat.location.start: feat.location.end - feat.shift]).reverse_complement())
                  else:
                    for x in feat.location.parts:
                      finSeq += str(rec.seq[feat.location.start + feat.shift: feat.location.end])
                  yield [
                    SeqRecord(
                        finSeq,
                        id=nid.replace(" ", "-"),
                        description=description,
                    )
                  ]
                elif feat.strand == -1:
                  yield [
                    SeqRecord(
                        (rec.seq[feat.location.start: feat.location.end - feat.shift]).reverse_complement(),
                        id=nid.replace(" ", "-"),
                        description=description,
                    )
                  ]
                else:
                  yield [
                    SeqRecord(
                        #feat.extract(rec).seq,
                        rec.seq[feat.location.start + feat.shift: feat.location.end],
                        id=nid.replace(" ", "-"),
                        description=description,
                    )
                  ]
            rec.features = newfeats
            rec.annotations = {}
            #gffWrite([rec], sys.stdout)
    else:
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        
        for rec in gffParse(gff3, base_dict=seq_dict):
            noMatch = True
            if "Alias" in rec.features[0].qualifiers.keys():
                lColumn = rec.features[0].qualifiers["Alias"][0]
            else:
                lColumn = ""
            for x in seq_dict:
                if x == rec.id or x == lColumn:
                    noMatch = False
            if noMatch:
                sys.stderr.write("Error: No Fasta ID matches GFF ID '" + rec.id + "'\n")
                exit(1)
            for feat in sorted(
                feature_lambda(
                    rec.features,
                    feature_test_type,
                    {"type": feature_filter},
                    subfeatures=True,
                ),
                key=lambda f: f.location.start,
            ):
                id = feat.id
                if len(id) == 0:
                    id = get_id(feat)

                if nodesc:
                    description = ""
                else:
                    if feat.strand == -1:
                      important_data = {"Location": FeatureLocation(feat.location.start + 1, feat.location.end - feat.shift, feat.strand)}
                    else:
                      important_data = {"Location": FeatureLocation(feat.location.start + 1 + feat.shift, feat.location.end, feat.strand)}
                    if "Name" in feat.qualifiers:
                        important_data["Name"] = feat.qualifiers.get("Name", [""])[0]

                    description = "[{}]".format(
                        ";".join(
                            [
                                "{key}={value}".format(key=k, value=v)
                                for (k, v) in important_data.items()
                            ]
                        )
                    )

                if isinstance(feat.location, CompoundLocation):
                  finSeq = ""
                  if feat.strand == -1:
                    for x in feat.location.parts:
                      finSeq += str((rec.seq[x.start: x.end - feat.shift]).reverse_complement())
                  else:
                    for x in feat.location.parts:
                      finSeq += str(rec.seq[x.start + feat.shift: x.end])
                  yield [
                    SeqRecord(
                        Seq(finSeq),
                        id=id.replace(" ", "-"),
                        description=description,
                    )
                  ]

                else:

                  if feat.strand == -1:
                    yield [
                      SeqRecord(
                          seq=Seq(str(rec.seq[feat.location.start: feat.location.end - feat.shift])).reverse_complement(),
                          id=id.replace(" ", "-"),
                          description=description,
                      )
                    ]
                  else:
                    yield [
                      SeqRecord(
                          #feat.extract(rec).seq,
                          seq=Seq(str(rec.seq[feat.location.start + feat.shift: feat.location.end])),
                          id=id.replace(" ", "-"),
                          description=description,
                      )
                    ]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Export corresponding sequence in genome from GFF3", epilog=""
    )
    parser.add_argument("fasta", type=argparse.FileType("r"), help="Fasta Genome")
    parser.add_argument("gff3", type=argparse.FileType("rb"), help="GFF3 File")
    parser.add_argument(
        "--feature_filter", default=None, help="Filter for specific feature types"
    )
    parser.add_argument(
        "--nodesc", action="store_true", help="Strip description field off"
    )
    args = parser.parse_args()
    for seq in main(**vars(args)):
        #if isinstance(seq, list):
        #  for x in seq:
        #    print(type(x.seq))
        #    SeqIO.write(x, sys.stdout, "fasta")
        #else:
          SeqIO.write(seq, sys.stdout, "fasta")
