#!/usr/bin/env python
import sys
import argparse
import logging
from intervaltree import IntervalTree
from Bio import SeqIO
from BCBio import GFF
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def extract_gff3_regions(gff3_files):
    data = {}
    for file in gff3_files:
        for record in GFF.parse(file):
            if record.id not in data:
                data[record.id] = IntervalTree()

            for gene in record.features:
                if gene.type in ('remark', 'DNA', 'annotation'):
                    continue

                data[record.id][int(gene.location.start):int(gene.location.end)] = True
    # Merge overlapping intervals
    for key in data:
        data[key].merge_overlaps()
    return data


def gaps(interval):
    sinterval = list(sorted(interval.items()))
    for i in range(len(sinterval) - 1):
        # Do we really care if we don't yield the wrap around one?
        yield (sinterval[i].end, sinterval[i + 1].begin)


def safe_reopen(fasta_file=None, gff3_files=None):
    occupied_regions = extract_gff3_regions(gff3_files)

    for record in SeqIO.parse(fasta_file, 'fasta'):
        # Get our list of gaps for this record
        gaps_in_data = list(gaps(occupied_regions[record.id]))
        # Arbitrarily choose the last one, so we re-open a bit upstream
        after = gaps_in_data[-1]
        # Midpoint
        after = sum(after) / 2
        record = record[after:] + record[0:after]
        record.description += ' [SafelyReopend=%s,%s bases upstream to avoid features]' % (after, after - len(record))
        yield record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta_file', type=argparse.FileType("r"))
    parser.add_argument('gff3_files', type=argparse.FileType("r"), nargs='+')

    args = parser.parse_args()
    for rec in safe_reopen(**vars(args)):
        SeqIO.write([rec], sys.stdout, "fasta")
