#!/usr/bin/env python
from Bio import SeqIO
import sys
from xmfa import parse_xmfa, percent_identity
import argparse
import logging
import itertools

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def _id_tn_dict(sequences):
    """Figure out sequence IDs AND sequence lengths from fasta file
    """
    label_convert = {}
    if sequences is not None:
        if len(sequences) == 1:
            for i, record in enumerate(SeqIO.parse(sequences[0], "fasta")):
                label_convert[str(i + 1)] = {"id": record.id, "len": len(record)}
        else:
            for i, sequence in enumerate(sequences):
                for record in SeqIO.parse(sequence, "fasta"):
                    label_convert[str(i + 1)] = {"id": record.id, "len": len(record)}
                    continue

    return label_convert


def total_similarity(xmfa_file, sequences=None, dice=False):
    if sequences is None:
        raise Exception("Must provide a non-zero number of sequence files")

    label_convert = _id_tn_dict(sequences)
    lcbs = parse_xmfa(xmfa_file)

    # make a matrix based on number of sequences
    table = {}

    for lcb in lcbs:
        # ignore LCBs containing only one sequence
        if len(lcb) == 0:
            continue

        # permutations based on num sequences to compare for current LCB
        compare_seqs = list(itertools.permutations(range(0, len(lcb)), 2))
        for permutation in compare_seqs:
            (i, j) = permutation
            similarity = percent_identity(lcb[i]["seq"], lcb[j]["seq"])

            i_name = label_convert[lcb[i]["id"]]["id"]
            j_name = label_convert[lcb[j]["id"]]["id"]
            # find length of sequence in LCB
            length_seq_lcb = lcb[i]["end"] - (lcb[i]["start"] - 1)
            # populate table with normalized similarity value based on length_seq_lcb
            if (i_name, j_name) not in table:
                table[(i_name, j_name)] = 0
            table[(i_name, j_name)] += length_seq_lcb * similarity

    # finalize total percent similarity by dividing by length of parent sequence
    for i in label_convert.keys():
        for j in label_convert.keys():
            i_name = label_convert[i]["id"]
            j_name = label_convert[j]["id"]
            if (i_name, j_name) in table:
                if dice:
                    table[(i_name, j_name)] = (
                        2
                        * table[(i_name, j_name)]
                        / (label_convert[i]["len"] + label_convert[j]["len"])
                    )
                else:
                    table[(i_name, j_name)] = (
                        table[(i_name, j_name)] / label_convert[i]["len"]
                    )
            else:
                table[(i_name, j_name)] = 0

            if i_name == j_name:
                table[(i_name, j_name)] = 100

    # print table
    names = []
    table_keys = sorted(label_convert.keys())

    for i in table_keys:
        names.append(label_convert[i]["id"])

    sys.stdout.write("\t" + "\t".join(names) + "\n")
    for j in table_keys:
        j_key = label_convert[j]["id"]
        sys.stdout.write(j_key)
        for i in table_keys:
            i_key = label_convert[i]["id"]
            sys.stdout.write("\t%0.2f" % table[(i_key, j_key)])
        sys.stdout.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert XMFA alignments to gff3", prog="xmfa2gff3"
    )
    parser.add_argument("xmfa_file", type=argparse.FileType("r"), help="XMFA File")
    parser.add_argument(
        "sequences",
        type=argparse.FileType("r"),
        nargs="+",
        help="Fasta files (in same order) passed to parent for reconstructing proper IDs",
    )
    parser.add_argument(
        "--dice", action="store_true", help="Use dice method for calculating % identity"
    )
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")

    args = parser.parse_args()

    total_similarity(**vars(args))
