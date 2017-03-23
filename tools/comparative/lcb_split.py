#!/usr/bin/env python
import argparse
import copy
import logging
import xmfa
from itertools import groupby
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def split_lcb(lcb, window_size=10, threshold=0.7):
    # Transpose sequence
    lines = []
    max_align_num = len(lcb[0]['seq'])
    for i in range(max_align_num):
        lines.append([])
        for j in range(len(lcb)):
            c = lcb[j]['seq'][i]
            if c != '-':
                lines[i].append(j)

    count_groups = []
    for i in range(0, len(lines), window_size):
        current_lines = lines[i:i+window_size]
        flat_list = [a for b in current_lines for a in b]
        counts = []
        for i in range(len(lcb)):
            value = float(flat_list.count(i)) / window_size
            if value >= threshold:
                counts.append(i)
        count_groups.append(counts)

    # groups = [(next(j), len(list(j)) + 1) for i, j in ]
    # [([4], 2), ([2, 3, 4, 5, 6], 2), ([0, 1, 2, 3, 4, 5, 6], 14), ([0, 3], 1)]
    # This says for 2 window sizes, we emit a new LCB with just [0:10] and
    # [10:20] for lcb #4, then one with all but 0/1 for 2, then all for 14.
    new_lcbs = []
    position = 0
    for i, j in groupby(count_groups):
        tmp = list(j)
        count = len(tmp)
        members = tmp[0]
        local_members = []
        for member in members:
            tmp_member = copy.deepcopy(lcb[member])
            tmp_member['seq'] = tmp_member['seq'][window_size * position:window_size * (position + count)]
            tmp_member['start'] = tmp_member['start'] + (window_size * position)
            tmp_member['end'] = tmp_member['start'] + (window_size * count)
            local_members.append(tmp_member)
        if len(local_members) > 0:
            new_lcbs.append(local_members)

        position += count
    return new_lcbs


def split_lcbs(lcbs, window_size=10, threshold=100):
    new_lcbs = []
    for lcb in lcbs:
        new_lcbs.extend(split_lcb(lcb, window_size=window_size, threshold=threshold))
    return new_lcbs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split XMFA alignments', prog='xmfa2smallerXmfa')
    parser.add_argument('xmfa_file', type=argparse.FileType("r"), help='XMFA File')

    parser.add_argument('--window_size', type=int, help='Window size for analysis', default=10)
    parser.add_argument('--threshold', type=float, help='All genomes must meet N percent similarity', default=0.7)

    args = parser.parse_args()

    # Write
    xmfa.to_xmfa(
        # Split
        split_lcbs(
            # Parse
            xmfa.parse_xmfa(args.xmfa_file),
            window_size=args.window_size,
            threshold=args.threshold,
        )
    )
