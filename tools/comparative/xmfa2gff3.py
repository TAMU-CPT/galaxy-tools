#!/usr/bin/env python
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
import argparse
from BCBio import GFF
import logging
from xmfa import parse_xmfa, _percent_identity
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
from xmfa import parse_xmfa, percent_identity, id_tn_dict


def generate_subfeatures(parent, window_size, other):
    log.debug("Generating subfeatures for %s %s %s", parent, window_size, other)
    for i in range(0, len(parent['seq']), window_size):
        block_seq = parent['seq'][i:i + window_size]
        real_window_size = len(block_seq)
        real_start = abs(parent['start']) - parent['seq'][0:i].count('-') + i - 1
        real_end = real_start + real_window_size - block_seq.count('-')

        log.debug("  I: %s, BS: %s, RWS: %s, RS: %s, RE: %s", i, block_seq, real_window_size, real_start, real_end)

        # if (real_end - real_start) < 10:
            # continue

        if parent['start'] < 0:
            strand = -1
        else:
            strand = 1

        pid = percent_identity(block_seq, other['seq'][i:i + real_window_size])
        # Ignore 0% identity sequences
        if pid == 0:
            continue

        yield SeqFeature(
            FeatureLocation(real_start, real_end),
            type="match_part", strand=strand,
            qualifiers={
                "source": "progressiveMauve",
                'score': pid
            }
        )


def convert_xmfa_to_gff3(xmfa_file, sequences=None, window_size=1000):
    label_convert = id_tn_dict(sequences)
    lcbs = parse_xmfa(xmfa_file)

    parent_records = {
        x: SeqRecord(Seq("ACTG", IUPAC.IUPACUnambiguousDNA), id=label_convert[x]['record_id'])
        for x in label_convert.keys()
    }

    for lcb_idx, lcb in enumerate(lcbs):
        import pprint; pprint.pprint(lcb)
        ids = [seq['id'] for seq in lcb]
        # Skip sequences that are JUST a single genome
        if len(ids) == 1:
            continue

        for parent_id in ids:
            parent = [seq for seq in lcb if seq['id'] == parent_id][0]
            others = [seq for seq in lcb if seq['id'] != parent_id]

            if parent['start'] == 0 and parent['end'] == 0:
                continue

            for o_idx, other in enumerate(others):
                # A feature representing a region of synteny between parent and the given other
                other_feature = SeqFeature(
                    FeatureLocation(parent['start'] - 1, parent['end']),
                    type="match", strand=parent['strand'],
                    qualifiers={
                        "source": "progressiveMauve",
                        "target": label_convert[other['id']]['record_id'],
                        "ID": 'm_%s_%s_%s' % (lcb_idx, o_idx, label_convert[other['id']]['record_id'])
                    }
                )

                for subfeature in generate_subfeatures(parent, window_size, other):
                    other_feature.sub_features.append(subfeature)

                parent_records[parent['id']].features.append(other_feature)

    for i in parent_records:
        yield [parent_records[i]]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert XMFA alignments to gff3', prog='xmfa2gff3')
    parser.add_argument('xmfa_file', type=file, help='XMFA File')
    parser.add_argument('--window_size', type=int, help='Window size for analysis', default=1000)
    parser.add_argument('--sequences', type=file,
                        help='Fasta files (in same order) passed to parent for reconstructing proper IDs')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    for result in convert_xmfa_to_gff3(**vars(args)):
        GFF.write(result, sys.stdout)
