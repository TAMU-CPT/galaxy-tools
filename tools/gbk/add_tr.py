#!/usr/bin/env python
import argparse
import sys
import copy
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def add_tr(genbank_file, end=None):
    record = list(SeqIO.parse(genbank_file, "genbank"))[0]
    cut_start = 0
    cut_end = 0

    # Find region of interest
    if end is None:
        repeat_regions = [x for x in record.features if x.type == 'repeat_region']
        if len(repeat_regions) == 1:
            cut_start = repeat_regions[0].location.start
            cut_end = repeat_regions[0].location.end
        elif len(repeat_regions) == 0:
            log.error("No repeats and end not specified")
            sys.exit(1)
        else:
            cut_start = repeat_regions[0].location.start
            cut_end = repeat_regions[0].location.end
            log.warning("Multiple repeat_regions found, using first")
    else:
        cut_start = 0
        cut_end = int(end)

    # Clone features in region of interest for duplication, this will actually
    # grab the repeat_region as well, conveniently
    clonefeats = [copy.deepcopy(x) for x in record.features if
                    cut_start <= x.location.start <= cut_end and
                    cut_start <= x.location.end <= cut_end]

    # For each cloned feature update the location
    for feat in clonefeats:
        feat.location += len(record.seq)
        if 'locus_tag' in feat.qualifiers:
            feat.qualifiers['locus_tag'][0] += '_rep'

        if 'note' in feat.qualifiers:
            feat.qualifiers['note'].append('CPT_TR dupcliated feature')
        else:
            feat.qualifiers['note'] = ['CPT_TR dupcliated feature']

    # Append extra sequence
    record.seq = record.seq + record.seq[cut_start:cut_end]
    # Add extended feature set
    record.features.extend(clonefeats)
    # Return
    return [record]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add Terminal Repeats to a genbank file')
    parser.add_argument('genbank_file', type=file, help='Genbank file')
    parser.add_argument('--end', type=str, help='End of Terminal Repeat Region. READ THE DIRECTIONS BELOW')

    args = parser.parse_args()
    SeqIO.write(add_tr(**vars(args)), sys.stdout, 'genbank')
