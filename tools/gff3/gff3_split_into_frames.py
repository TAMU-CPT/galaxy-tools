#!/usr/bin/env python
import argparse
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
import copy


def split_into_frames(gff3):
    for rec in GFF.parse(gff3):
        rf1 = []
        rf2 = []
        rf3 = []
        rf4 = []
        rf5 = []
        rf6 = []
        dummy_rec = copy.deepcopy(rec)
        dummy_rec.annotations = {}
        for gene in feature_lambda(
            rec.features,
            feature_test_type,
            {'types': 'gene'},
            subfeatures=True
        ):
            if gene.location.strand == 1:
                frame = str(((gene.location.start) % 3) + 1)
            else:
                frame = str((-(gene.location.start - 1) % 3) + 4)
            locals()['rf' + frame].append(gene)

        for i in range(6):
            dummy_rec.features = locals()['rf' + str(i + 1)]
            with open('rf' + str(i + 1) + '.gff3', 'a') as outfile:
                GFF.write([dummy_rec], outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='split genes into separate files based on reading frame')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='GFF3 annotations')
    args = parser.parse_args()
    split_into_frames(**vars(args))
