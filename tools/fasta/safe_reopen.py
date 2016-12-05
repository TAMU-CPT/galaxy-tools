#!/usr/bin/env python
import sys
import argparse
import logging
from Bio import SeqIO
from BCBio import GFF
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def extract_gff3_regions(gff3_files):
    data = {}
    for idx, file in enumerate(gff3_files):
        for record in GFF.parse(file):
            if record.id not in data:
                data[record.id] = {}

            if idx not in data[record.id]:
                data[record.id][idx] = []

            sorted_genes = sorted(record.features, key=lambda feature: feature.location.start)
            if len(sorted_genes) == 0:
                continue

            current_gene = [
                int(sorted_genes[0].location.start),
                int(sorted_genes[0].location.end)
            ]
            for gene in sorted_genes[1:]:
                # If the gene's start is contiguous to the "current_gene", then we
                # extend current_gene
                if gene.location.start <= current_gene[1] + 10:
                    current_gene[1] = int(gene.location.end)
                else:
                    # If it's discontiguous, we append the region and clear.
                    data[record.id][idx].append(current_gene)
                    current_gene = [int(gene.location.start), int(gene.location.end)]

            # This generally expected that annotations would NOT continue unto the end
            # of the genome, however that's a bug, and we can make it here with an
            # empty contiguous_regions list
            data[record.id][idx].append(current_gene)
    return data


def gap_in_region(query, track):
    for region in track:
        if region[0] <= query <= region[1]:
            return False
    return True


def gap_in_data(data, record_length):
    # We'll pick the 0th track just for our reference.
    # We'll also skip the first couple because...yeah.
    for i in range(2, len(data[0]) + 1):
        if i == 0:
            a = (1, 1)
            b = data[0][i]
        elif i >= len(data[0]):
            a = data[0][i - 1]
            b = (record_length, None)
        else:
            a = data[0][i - 1]
            b = data[0][i]

        if abs(b[0] - a[1]) > 30:
            mid = float(b[0] + a[1]) / 2
            if len(data.keys()) == 1:
                return mid
            else:
                if all([
                    gap_in_region(mid, data[j]) for j in range(1, len(data.keys()))
                ]):
                    return int(mid)
    raise Exception("Could not find acceptable gap")


def safe_reopen(fasta_file=None, gff3_files=None):
    gff3_pos = extract_gff3_regions(gff3_files)
    for record in SeqIO.parse(fasta_file, 'fasta'):
        after = gap_in_data(gff3_pos[record.id], len(record))
        seq = record.seq
        record.seq = seq[after:] + seq[0:after]
        yield record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta_file', type=argparse.FileType("r"))
    parser.add_argument('gff3_files', type=argparse.FileType("r"), nargs='+')

    args = parser.parse_args()
    for rec in safe_reopen(**vars(args)):
        SeqIO.write([rec], sys.stdout, "fasta")
