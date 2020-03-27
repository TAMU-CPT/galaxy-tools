#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
from Bio import SeqIO
import re
import argparse
from phantasm import Cassettes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Cassette IDs for a set of genomes"
    )
    parser.add_argument(
        "genbank_file", type=argparse.FileType("r"), help="GenBank File"
    )
    parser.add_argument("--version", action="version", version="0.2")
    args = parser.parse_args()
    cassette_model = Cassettes.get_clustering_data()
    regex_containers = {}
    custom_regexes = {}

    for key in cassette_model:
        for member in cassette_model[key]["members"]:
            compiled = re.compile(member, re.IGNORECASE)
            regex_containers[member] = {
                "re": compiled,
                "parent": cassette_model[key]["id"],
            }
        if "custom" in cassette_model[key]:
            for custom in cassette_model[key]["custom"]:
                custom_regexes[custom] = {
                    "title": cassette_model[key]["id"],
                    "is": [
                        re.compile(x, re.IGNORECASE)
                        for x in cassette_model[key]["custom"][custom]["is"]
                    ],
                    "isnot": [
                        re.compile(x, re.IGNORECASE)
                        for x in cassette_model[key]["custom"][custom]["isnot"]
                    ],
                }

    print "\t".join(["#ID", "CID"])
    for i, record in enumerate(SeqIO.parse(args.genbank_file, "genbank")):
        data = ""
        for feature in [
            x for x in record.features if x.type == "CDS" and "product" in x.qualifiers
        ]:
            matched_result = None
            feature_result = "".join(feature.qualifiers["product"])

            for regex in regex_containers:
                if regex_containers[regex]["re"].match(feature_result):
                    matched_result = regex_containers[regex]["parent"]

            for custom in custom_regexes:
                care = False
                # If we hit to an is-statement, we care
                for regex in custom_regexes[custom]["is"]:
                    if regex.match(feature_result):
                        care = True
                # As long as we don't hit to an isnot-statement
                for regex in custom_regexes[custom]["isnot"]:
                    if regex.match(feature_result):
                        care = False

                if care:
                    is_okay = True
                    ok_to_overwrite = True

                    for regex in custom_regexes[custom]["is"]:
                        if not regex.match(feature_result):
                            is_okay = False

                    for regex in custom_regexes[custom]["isnot"]:
                        if regex.match(feature_result):
                            ok_to_overwrite = False
                            is_okay = False

                    if is_okay and ok_to_overwrite:
                        matched_result = custom_regexes[custom]["title"]

            if matched_result is not None:
                if feature.location.strand > 0:
                    data += "+"
                else:
                    data += "-"
                data += matched_result

        if len(data) > 4:
            mod = Cassettes.collapse(data)
            if mod is not None:
                print "\t".join([record.id, mod])
