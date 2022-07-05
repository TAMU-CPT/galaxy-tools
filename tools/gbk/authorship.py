#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
import logging

logging.basicConfig(level=logging.INFO)
REF_KEYS = ["authors", "title", "journal", "medline", "pubmed", "comment"]


def lineage(genbank_files=None, first=False):
    from Bio import SeqIO

    for genbank_file in genbank_files:
        records = list(SeqIO.parse(genbank_file, "genbank"))
        for record in records:
            id = record.id
            if "." in id:
                id = id.split(".")[0]

            if "references" in record.annotations:
                refs = record.annotations["references"]
                if first:
                    yield [id] + [
                        "" if not hasattr(refs[0], x) else getattr(refs[0], x)
                        for x in REF_KEYS
                    ]
                else:
                    for ref in refs:
                        yield [id] + [
                            "" if not hasattr(ref, x) else getattr(ref, x)
                            for x in REF_KEYS
                        ]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract authorship " "information from genbank file"
    )
    parser.add_argument(
        "genbank_files", type=argparse.FileType("r"), nargs="+", help="Genbank file"
    )
    parser.add_argument(
        "--first", action="store_true", help="Only pick out first reference"
    )

    parser.add_argument("--version", action="version", version="0.1")
    args = parser.parse_args()
    print ("\t".join(["# sequence", "reference"] + REF_KEYS))
    for line in lineage(**vars(args)):
        print ("\t".join(line))
