#!/usr/bin/env python
import os
import copy
import argparse
from BCBio import GFF
import logging

logging.basicConfig(level=logging.INFO)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=argparse.FileType("r"), help="GFF3 File")
    parser.add_argument("keys", type=str, nargs="+", help="Unique properties.")
    parser.add_argument(
        "--joiner", help="String to use in joining properties for output file names."
    )

    args = parser.parse_args()

    file_handles = {}
    logging.info("Parsing Data")
    for record in GFF.parse(args.data):
        logging.info("Process record %s", record.id)
        record_features = record.features
        record.features = []
        tmprec = copy.deepcopy(record)
        tmprec.annotations = {}
        tmprec.features = []
        record.features = record_features

        for feature in record.features:
            props = []
            if "record_id" in args.keys:
                props.append(record.id)
            if "source" in args.keys:
                props.append(feature.qualifiers["source"][0])
            if "target" in args.keys:
                props.append(feature.qualifiers["Target"][0])

            propkey = "|".join(map(str, props))

            if propkey not in file_handles:
                filename = args.joiner.join(props)
                path = os.path.join("out", filename + ".gff3")
                logging.info("Opening %s", path)
                file_handles[propkey] = open(path, "a")

            tmprec.features = [feature]
            GFF.write([tmprec], file_handles[propkey])
        # SeqIO.write([record], args.fasta, 'fasta')
        # sys.exit()

    for key in file_handles:
        file_handles[key].close()
