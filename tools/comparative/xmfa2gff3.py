#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse
from BCBio import GFF
import logging
from xmfa import parse_xmfa, _percent_identity
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def _id_tn_dict(sequences):
    """Figure out sequence IDs
    """
    label_convert = {}
    if sequences is not None:
        if len(sequences) == 1:
            for i, record in enumerate(SeqIO.parse(sequences[0], 'fasta')):
                label_convert[str(i + 1)] = record.id
        else:
            for i, sequence in enumerate(sequences):
                for record in SeqIO.parse(sequence, 'fasta'):
                    label_convert[str(i + 1)] = record.id
                    continue
    return label_convert


def convert_xmfa_to_gff3(xmfa_file, relative_to='1', sequences=None, window_size=1000):
    label_convert = _id_tn_dict(sequences)

    lcbs = parse_xmfa(xmfa_file)

    for lcb_idx, lcb in enumerate(lcbs):
        record = SeqRecord(Seq("A"), id=label_convert.get(relative_to, relative_to))
        ids = [seq['id'] for seq in lcb]

        # Doesn't match part of our sequence
        if relative_to not in ids:
            continue

        # Skip sequences that are JUST our "relative_to" genome
        if len(ids) == 1:
            continue

        parent = [seq for seq in lcb if seq['id'] == relative_to][0]
        others = [seq for seq in lcb if seq['id'] != relative_to]

        if parent['start'] == 0 and parent['end'] == 0:
            continue

        for o_idx, other in enumerate(others):
            other['feature'] = SeqFeature(
                FeatureLocation(parent['start'], parent['end'] + 1),
                type="match", strand=parent['strand'],
                qualifiers={
                    "source": "progressiveMauve",
                    "target": label_convert.get(other['id'], other['id']),
                    "ID": 'm_%s_%s_%s' % (lcb_idx, o_idx, label_convert.get(other['id'], 'xmfa_' + other['rid']))
                }
            )

        for i in range(0, len(lcb[0]['seq']), window_size):
            block_seq = parent['seq'][i:i + window_size]
            real_window_size = len(block_seq)
            real_start = abs(parent['start']) - parent['seq'][0:i].count('-') + i
            real_end = real_start + real_window_size - block_seq.count('-')

            if (real_end - real_start) < 10:
                continue

            if parent['start'] < 0:
                strand = -1
            else:
                strand = 1

            for other in others:
                pid = _percent_identity(block_seq, other['seq'][i:i + real_window_size])
                # Ignore 0% identity sequences
                if pid == 0:
                    continue
                other['feature'].sub_features.append(
                    SeqFeature(
                        FeatureLocation(real_start, real_end),
                        type="match_part", strand=strand,
                        qualifiers={
                            "source": "progressiveMauve",
                            'score': pid
                        }
                    )
                )

        for other in others:
            record.features.append(other['feature'])
        record.annotations = {}
        yield [record]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert XMFA alignments to gff3', prog='xmfa2gff3')
    parser.add_argument('xmfa_file', type=file, help='XMFA File')
    parser.add_argument('--window_size', type=int, help='Window size for analysis', default=1000)
    parser.add_argument('--relative_to', type=str, help='Index of the parent sequence in the MSA', default='1')
    parser.add_argument('--sequences', type=file, nargs='+',
                        help='Fasta files (in same order) passed to parent for reconstructing proper IDs')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    for result in convert_xmfa_to_gff3(**vars(args)):
        GFF.write(result, sys.stdout)
