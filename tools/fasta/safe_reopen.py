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

            # Mark off the end of the genome as not available.
            data[record.id][len(record):(len(record) + 1)] = True
    # Merge overlapping intervals
    for key in data:
        data[key].merge_overlaps()
    return data


def gaps(interval):
    sinterval = list(sorted(interval.items()))
    for i in range(len(sinterval) - 1):
        # Do we really care if we don't yield the wrap around one?
        yield (sinterval[i].end, sinterval[i + 1].begin)


def nearest_gap(gaps, position, strand):
    # TODO: Cannot handle starting at base 1 and looking upstream
    if position == -2:
        after = gaps[-1]
        return sum(after) / 2
    elif position == -1:
        after = gaps[0]
        return sum(after) / 2
    else:
        best = None
        if strand > 0:
            for gap in sorted(gaps, key=lambda x: x[0]):
                if gap[1] < position:
                    if not best:
                        best = gap

                    if position - gap[1] < position - best[1]:
                        best = gap
                elif gap[0] < position < gap[1]:
                    best = gap
                    break
        else:
            for gap in sorted(gaps, key=lambda x: -x[0]):
                if gap[0] > position:
                    if not best:
                        best = gap

                    if gap[0] - position < best[0] - position:
                        best = gap

                elif gap[0] < position < gap[1]:
                    best = gap
                    break
        return sum(best) / 2


def safe_reopen(fasta_file=None, gff3_files=None, position=-1, strand=0):
    occupied_regions = extract_gff3_regions(gff3_files)

    for record in SeqIO.parse(fasta_file, 'fasta'):
        # Get our list of gaps for this record
        gaps_in_data = list(gaps(occupied_regions[record.id]))
        # Arbitrarily choose the last one, so we re-open a bit upstream
        after = nearest_gap(gaps_in_data, position, strand)
        # Midpoint
        record = record[after:] + record[0:after]
        # If it's a minus strand, auto-revcom
        if strand == -1:
            record.seq = record.seq.reverse_complement()
            record.description += ' [SafeReopenRevCom=True]'

        record.description += ' [SafelyReopened=%s,%s bases from end]' % (after, len(record) - after)
        yield record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta_file', type=argparse.FileType("r"))
    parser.add_argument('gff3_files', type=argparse.FileType("r"), nargs='+')
    parser.add_argument('--position', type=int, default=-1)
    parser.add_argument('--strand', type=int, default=1)

    args = parser.parse_args()
    for rec in safe_reopen(**vars(args)):
        SeqIO.write([rec], sys.stdout, "fasta")
