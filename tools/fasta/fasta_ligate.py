#!/usr/bin/env python
import sys
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def merge(sequences, id="merged", pre="", post="", link="", bed="out.bed"):
    seqs = []
    ids = []
    for sequence in sequences:
        for record in SeqIO.parse(sequence, "fasta"):
            seqs.append(str(record.seq))
            ids.append(record.id)

    with open(bed, "w") as bed_output:
        start_offset = len(pre)
        bed_tpl = "{chr}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"
        for (s, i) in zip(seqs, ids):
            bed_output.write(
                bed_tpl.format(
                    chr=id,
                    start=start_offset,
                    end=start_offset + len(s),
                    name=i,
                    score=1000,
                    strand="+",
                )
            )
            start_offset += len(s) + len(link)

    output_seq = pre + link.join(seqs) + post

    output = SeqRecord(
        seq=Seq(output_seq),
        id=id,
        description="Created from [%s...]" % (",".join(ids[0:10])),
    )
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ligate sequences")
    parser.add_argument(
        "--id", help="New fasta identifier for merged sequences", default="merged"
    )

    parser.add_argument("--pre", help="Pre-sequence addition")
    parser.add_argument("--post", help="Post-sequence addition")
    parser.add_argument("--link", help="Linker sequence")
    parser.add_argument("--bed", help="Output bed data")

    parser.add_argument(
        "sequences",
        type=argparse.FileType("r"),
        nargs="+",
        help="Fasta sequences (1 or more per file)",
    )

    args = parser.parse_args()
    SeqIO.write(merge(**vars(args)), sys.stdout, "fasta")
