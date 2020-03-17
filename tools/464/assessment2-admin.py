#!/usr/bin/env python
import logging
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
from guanine import GuanineClient

logging.basicConfig(level=logging.INFO)
logging.getLogger("requests").setLevel(logging.CRITICAL)
log = logging.getLogger(__name__)
g = GuanineClient()


def validate(ogs, user_gff3, user_email, offset=213):
    comp = {}
    for rec in GFF.parse(ogs):
        for feature in feature_lambda(
            rec.features, feature_test_type, {"type": "CDS"}, subfeatures=True
        ):
            if feature.strand > 0:
                offset_end = int(feature.location.end)
            else:
                offset_end = int(feature.location.start)

            comp[offset_end] = feature
    max_score = len(comp)

    user = {}
    for rec in GFF.parse(user_gff3):
        for feature in feature_lambda(
            rec.features, feature_test_type, {"type": "gene"}, subfeatures=True
        ):
            cdss = list(
                feature_lambda(
                    feature.sub_features,
                    feature_test_type,
                    {"type": "CDS"},
                    subfeatures=True,
                )
            )
            if len(cdss) == 0:
                continue

            cds = cdss[0]
            cds.qualifiers["Name"] = feature.qualifiers.get("Name", [])

            if cds.strand > 0:
                offset_end = int(cds.location.end) + offset
            else:
                offset_end = int(cds.location.start) + offset
            user[offset_end] = cds

    results = []

    for user_annotation in sorted(user.keys()):
        fid = user[user_annotation].qualifiers.get("Name", [user[user_annotation].id])[
            0
        ]
        if user_annotation in comp:
            # User successfully annotated gene
            # print comp[user_annotation]
            # good.append(user_annotation)
            ogs_feature = comp[user_annotation]
            usr_feature = user[user_annotation]
            if (
                ogs_feature.location.start == usr_feature.location.start + offset
                and ogs_feature.location.end == usr_feature.location.end + offset
            ):
                # Success!
                results.append(
                    {
                        "points": 1,
                        "message": "Correct",
                        "stop": user_annotation - offset,
                        "id": fid,
                    }
                )
            else:
                results.append(
                    {
                        "points": 0.5,
                        "message": "Wrong start codon",
                        "stop": user_annotation - offset,
                        "id": fid,
                    }
                )
            del comp[user_annotation]
        else:
            # User annotated a non-existent gene
            results.append(
                {
                    "points": 0,
                    "message": "Not an actual gene, according to the official gene set",
                    "stop": user_annotation - offset,
                    "id": fid,
                }
            )

    for leftover in comp.keys():
        results.append(
            {
                "points": -0.5,
                "message": "Missed a gene in the official gene set",
                "stop": leftover,
                "id": "---",
            }
        )

    score = sum(x["points"] for x in results)
    print("Final score: %s / %s" % (score, max_score))
    print()

    # Submit score to GUANINE
    final_score = float(float(score) / max_score)
    if final_score < 0:
        final_score = 0
    g.submit(user_email, "C1", final_score)

    for x in results:
        print(
            "Feature {id} ending with stop codon: {stop}; Points: {points}; Message {message}".format(
                **x
            )
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="verify against expectations")
    parser.add_argument("ogs", type=argparse.FileType("r"), help="GFF3 annotations")
    parser.add_argument(
        "user_gff3", type=argparse.FileType("r"), help="GFF3 annotations"
    )
    parser.add_argument("user_email", help="User email")
    args = parser.parse_args()
    validate(**vars(args))
