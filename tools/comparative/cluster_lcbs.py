#!/usr/bin/env python
from Bio import SeqIO
import tempfile
import sys
import argparse
import operator


def parse_xmfa(xmfa):
    """Simple XMFA parser until https://github.com/biopython/biopython/pull/544
    """
    current_lcb = []
    current_seq = {}
    for line in xmfa.readlines():
        if line.startswith('#'):
            continue

        if line.strip() == '=':
            if 'id' in current_seq:
                current_lcb.append(current_seq)
                current_seq = {}
            yield current_lcb
            current_lcb = []
        else:
            line = line.strip()
            if line.startswith('>'):
                if 'id' in current_seq:
                    current_lcb.append(current_seq)
                    current_seq = {}
                data = line.strip().split()
                # 0 1           2 3      4 5
                # > 1:5986-6406 + CbK.fa # CbK_gp011
                id, loc = data[1].split(':')
                start, end = loc.split('-')
                current_seq = {
                    'rid': '_'.join(data[1:]),
                    'id': id,
                    'start': int(start),
                    'end': int(end),
                    'strand': 1 if data[2] == '+' else -1,
                    'file': data[3],
                    'seq': '',
                    'comment': '',
                }
                if len(data) > 5:
                    current_seq['comment'] = ' '.join(data[5:])
            else:
                current_seq['seq'] += line.strip()


HEADER_TPL = '> {id}:{start}-{end} {strand} {file} # {comment}\n'


def split_by_n(seq, n):
    """A generator to divide a sequence into chunks of n units."""
    # http://stackoverflow.com/questions/9475241/split-python-string-every-nth-character
    while seq:
        yield seq[:n]
        seq = seq[n:]


def to_xmfa(lcbs, handle=sys.stdout):
    handle.write("#FormatVersion Mauve1\n")
    for lcb in lcbs:
        for aln in lcb:
            handle.write(HEADER_TPL.format(
                id=aln['id'],
                start=aln['start'],
                end=aln['end'],
                strand='+' if aln['strand'] > 0 else '-',
                file=aln['file'],
                comment=aln['comment'],
            ))

            for line in split_by_n(aln['seq'], 80):
                handle.write(line + '\n')
        handle.write('=\n')


def percent_identity(a, b):
    """Calculate % identity, ignoring gaps in the host sequence
    """
    match = 0
    mismatch = 0
    for char_a, char_b in zip(list(a), list(b)):
        if char_a == '-':
            continue
        if char_a == char_b:
            match += 1
        else:
            mismatch += 1

    if match + mismatch == 0:
        return 0.0
    return 100 * float(match) / (match + mismatch)


def id_tn_dict(sequences, tmpfile=False):
    """Figure out sequence IDs
    """
    label_convert = {}
    correct_chrom = None
    if not isinstance(sequences, list):
        sequences = [sequences]

    i = 0
    for sequence_file in sequences:
        for record in SeqIO.parse(sequence_file, 'fasta'):
            if correct_chrom is None:
                correct_chrom = record.id

            i += 1
            key = str(i)
            label_convert[key] = {
                'record_id': record.id,
                'len': len(record.seq),
            }

            if tmpfile:
                label_convert[key] = tempfile.NamedTemporaryFile(delete=False)

    return label_convert

def filter_lcbs_for_seq(xmfa):
    """ clusters lcbs based on which sequences they involve """
    clusters = {}

    for i in list(parse_xmfa(xmfa)):
        cluster_name = ''

        for g in i:
            cluster_name += g['id']

        if cluster_name not in clusters:
            clusters[cluster_name] = [i]
        else:
            clusters[cluster_name].append(i)

    return clusters

    # to_xmfa(clusters['123456'])

def merge_lcbs(lcb1, lcb2):
    for num, i in enumerate(lcb1):
        i['end'] = lcb2[num]['end']
        i['seq'] += lcb2[num]['seq']

    return lcb1

def resolve_clusters(clusters):
    merged = []
    for lcbs in clusters:
        if len(lcbs) == 1:
            merged.append(lcbs[0])
            continue
        merging = lcbs[0]
        for lcb in lcbs[1:]:
            merging = merge_lcbs(merging, lcb)
        merged.append(merging)

    return merged

def new(clusters, lcb):
    new = True
    for c in clusters:
        if lcb in c:
            new = False
    return new

def cluster_lcbs(lcbs, threshold):
    """ clusters lcbs based on how far apart they are"""

    clusters = []
    for o, i in enumerate(lcbs):
        cluster = []

        if not new(clusters, i):
            continue

        cluster.append(i)
        compare_against = i

        for n, j in enumerate(lcbs):

            if not new(clusters, j) or i == j or compare_against == j:
                continue

            close = True
            for num, k in enumerate(compare_against):
            # for num, k in enumerate(i):
                if j[num]['start'] - k['end'] > threshold:
                    close = False

            if close:
                cluster.append(j)
                compare_against = j

        clusters.append(cluster)

    return resolve_clusters(clusters)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='process XMFA')
    parser.add_argument('xmfa', type=argparse.FileType("r"), help='XMFA file')
    parser.add_argument('threshold', type=int, help='maximum number of nucleotides between lcbs in a cluster')
    args = parser.parse_args()

    # assuming lcbs are filtered
    final_lcbs = []
    lcbs_filtered_for_seq = filter_lcbs_for_seq(args.xmfa)
    for i in lcbs_filtered_for_seq:
        final_lcbs += cluster_lcbs(lcbs_filtered_for_seq[i], args.threshold)
    to_xmfa(final_lcbs)
