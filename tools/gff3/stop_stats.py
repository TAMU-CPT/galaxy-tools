#!/usr/bin/env python
import argparse
from cpt_gffParser import gffParse, gffWrite
from Bio import SeqIO
from gff3 import feature_lambda, feature_test_type


def main(fasta, gff3):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    codon_usage = {}

    for rec in gffParse(gff3, base_dict=seq_dict):
        for feat in feature_lambda(
            rec.features, feature_test_type, {"type": "CDS"}, subfeatures=True
        ):
            seq = str(feat.extract(rec).seq)[-3:]
            try:
                codon_usage[seq] += 1
            except KeyError:
                codon_usage[seq] = 1

    names = {
        "TAG": "Amber",
        "TAA": "Ochre",
        "TGA": "Opal",
    }

    # TODO: print all actg combinations? Or just ones that are there
    print "# Name\tCodon\tCount"
    for key in sorted(codon_usage):
        print "\t".join((names.get(key.upper(), "None"), key, str(codon_usage[key])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Summarise stop codon usage", epilog=""
    )
    parser.add_argument("fasta", type=argparse.FileType("r"), help="Fasta Genome")
    parser.add_argument("gff3", help="GFF3 File")
    args = parser.parse_args()
    main(**vars(args))
