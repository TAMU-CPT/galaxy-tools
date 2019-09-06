#!/usr/bin/env python
import sys
import logging
import argparse
from Bio import SeqIO
from Bio.Data import CodonTable

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def translate(fasta_file, target="protein", table=11, strip_stops=False):
    records = list(SeqIO.parse(fasta_file, "fasta"))

    for record in records:
        if target == "protein":
            mod = len(record.seq) % 3
            if mod != 0:
                record.seq = record.seq[0:-mod]

            # Read http://biopython.org/DIST/docs/api/Bio.Seq.Seq-class.html#transcribe
            # for valid CDS conditions.

            # Will first try to translate sequence as a CDS,
            # then just as a sequence if this fails.

            try:
                tmpseq = record.seq.translate(table=table, cds=True)
            except CodonTable.TranslationError as cte:
                log.info("Translation issue at %s: %s", record.id, cte)
                tmpseq = record.seq.translate(table=table, cds=False)

            # check if stop in middle of protein
            if "*" in tmpseq:
                log.info(
                    "Trimming %s from %s to %s due to stop codons",
                    record.id,
                    len(record.seq),
                    3 * len(tmpseq) - 3,
                )
                tmpseq = tmpseq[0 : str(tmpseq).index("*")]

            # add stop to end if strip_stops=False
            if not strip_stops:
                tmpseq = tmpseq + "*"

            record.seq = tmpseq
            if len(record.seq) > 0:
                SeqIO.write(record, sys.stdout, "fasta")
        else:
            record.seq = record.seq.transcribe()
            SeqIO.write(record, sys.stdout, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate fasta file")
    parser.add_argument("fasta_file", type=argparse.FileType("r"), help="Fasta file")
    parser.add_argument("--target", choices=["protein", "rna"])
    parser.add_argument(
        "--table",
        type=int,
        default=11,
        help="Translation table to use",
        choices=range(1, 23),
    )
    parser.add_argument(
        "--strip_stops", action="store_true", help="Remove stop characters"
    )

    args = parser.parse_args()
    translate(**vars(args))
