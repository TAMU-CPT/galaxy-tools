#!/usr/bin/env python
import re
import BIO_FIX_TOPO  # NOQA
import os
import argparse
from Bio import SeqIO

import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

PHAGE_IN_MIDDLE = re.compile("^(?P<host>.*)\s*(phage|virus) (?P<phage>.*)$")
BACTERIOPHAGE_IN_MIDDLE = re.compile("^(?P<host>.*)\s*bacteriophage (?P<phage>.*)$")
STARTS_WITH_PHAGE = re.compile(
    "^(bacterio|vibrio|Bacterio|Vibrio)?[Pp]hage (?P<phage>.*)$"
)
NEW_STYLE_NAMES = re.compile("^(?P<phage>v[A-Z]_[A-Z][a-z]{2}[A-Z]_.*)$")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split GBK file")
    parser.add_argument(
        "genbank_files", type=argparse.FileType("r"), nargs="+", help="Genbank files"
    )
    parser.add_argument(
        "--style",
        type=str,
        choices=["host-first", "phage-first"],
        default="host-first",
        help="Set style format of new name. host-first or phage-first.",
    )
    args = parser.parse_args()

    outdir = os.path.join(os.getcwd(), "gbk_out")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for genbank_file in args.genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            pn = record.annotations["source"]
            m = PHAGE_IN_MIDDLE.match(pn)
            host = None
            phage = None
            if m:
                host = m.group("host")
                phage = m.group("phage")
            m = BACTERIOPHAGE_IN_MIDDLE.match(pn)
            if m:
                host = m.group("host")
                phage = m.group("phage")

            m = STARTS_WITH_PHAGE.match(pn)
            if m:
                phage = m.group("phage")

            m = NEW_STYLE_NAMES.match(pn)
            if m:
                phage = m.group("phage")

            if host and phage:
                host = host.strip()
                phage = phage.strip()
                if args.style == "host-first":
                    name = os.path.join(outdir, "%s phage %s.gbk" % (host, phage))
                else:
                    name = os.path.join(outdir, "%s %s phage.gbk" % (phage, host))
            elif phage:
                phage = phage.strip()
                if args.style == "host-first":
                    name = os.path.join(outdir, "Phage %s.gbk" % phage)
                else:
                    name = os.path.join(outdir, "%s.gbk" % phage)
            else:
                log.info("Could not determine name from %s. Contact IT.", pn)
                name = os.path.join(outdir, "%s.gbk" % record.id)

            with open(name, "w") as handle:
                log.info("Storing %s to %s", record.id, name)
                SeqIO.write([record], handle, "genbank")
