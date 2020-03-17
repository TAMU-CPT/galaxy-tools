#!/usr/bin/env python
import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RevCom fasta files")
    parser.add_argument(
        "fa_files", type=argparse.FileType("r"), nargs="+", help="fasta files"
    )

    kwargs = dict(
        id=True,
        name=True,
        description=True,
        features=True,
        annotations=True,
        letter_annotations=True,
        dbxrefs=True,
    )

    args = parser.parse_args()

    i = 1

    for fa_file in args.fa_files:
        # print(dir(fa_file))
        if fa_file.readline().startswith(">"):  # Determine if fasta
            fa_file.seek(0)  # Reset input buffer to start (Undo our readline)
            for record in SeqIO.parse(fa_file, "fasta"):
                SeqIO.write(record.reverse_complement(**kwargs), sys.stdout, "fasta")
        else:
            fa_file.seek(0)
            name = "Txt_File_Input_" + str(i)
            i += 1
            tempSeq = Seq("")
            for line in fa_file.readlines():
                tempSeq = tempSeq + line.strip("\n")  # .add(Seq(line.strip('\n')))
            tempRec = SeqRecord(
                tempSeq, name, name, "Reverse complement from raw .txt input."
            )
            SeqIO.write(tempRec.reverse_complement(**kwargs), sys.stdout, "fasta")
