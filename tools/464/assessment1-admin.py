#!/usr/bin/env python
import argparse
import json
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
from guanine_report import auth, post_result, student_id

# from guanine import GuanineClient
GUANINE_URL = "https://cpt.tamu.edu/guanine-backend/"


def validate(gff3):
    results = {}
    for rec in GFF.parse(gff3):
        for feature in feature_lambda(
            rec.features, feature_test_type, {"type": "gene"}, subfeatures=True
        ):
            checks = []
            graded = []
            # dbxrefs
            if "CPT:283675" in feature.qualifiers.get("Dbxref", []):
                checks.append(True)
                graded.append({})
            else:
                checks.append(False)
                graded.append({"q1": "0"})  # ???

            # Notes
            if "Howdy!" in feature.qualifiers.get("Note", []):
                checks.append(True)
                graded.append({})
            else:
                checks.append(False)
                graded.append({"q2": "0"})  # ???

            owner = feature.qualifiers.get("owner", ["unknown"])[0]
            results[owner] = {
                "checks": checks,
                "graded": graded,
                "score": checks.count(True),
            }

    # Process all students at once
    token = auth(open("/galaxy/creds.json", "r"), GUANINE_URL)
    for email, result in results.items():
        sid = student_id(email, GUANINE_URL, token)
        result = post_result(
            sid,
            result["score"],
            2,
            token,
            GUANINE_URL,
            "a59a5001-57e7-4776-8807-63b544735f3f",
            json.dumps({"raw": result, "graded": result["graded"]}),
        )
        if result.status_code in (200, 201):
            print("Success")
        else:
            print("[Error] user=%s msg=%s" % (email, result.text))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="verify against expectations")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    validate(**vars(args))
