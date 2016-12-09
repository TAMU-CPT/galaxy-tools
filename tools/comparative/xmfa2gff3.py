#!/usr/bin/env python
import os
import sys
import argparse
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from BCBio import GFF
from xmfa import parse_xmfa, percent_identity, id_tn_dict
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def generate_subfeatures(parent, window_size, other, protein=False):
    log.debug("Generating subfeatures for %s %s %s", parent, window_size, other)
    # We strip off trailing mismatches, and add a stop codon.
    #
    # The trailing mismatches cause blocks to run past the end where we don't
    # want them (and violate the gff3 spec).
    #
    # The stop codon will never match (and will slightly decrease %ident) but
    # it makes the last block the same length as the others which looks nicer.
    ps = parent['seq'].rstrip('-') + '*'
    for i in range(0, len(ps), window_size):
        block_seq = ps[i:i + window_size]
        real_window_size = len(block_seq)
        real_start = abs(parent['start'])
        # This is the index against the parent feature start/end region where our current comparison block is
        #
        # In order to handle cases where the parent sequence has a bunch of
        # ----- leading, we want to find the "real" start of the current block.
        # So we take the current block start (i), and subtract the number of -s
        #
        # Then convert to nucl.
        nucl2prot = 3
        indexed_start = nucl2prot * (i - ps[0:i].count('-'))

        log.debug("  I: %s, BS: %s, RWS: %s, RS: %s, IS: %s, C: %s", i,
                  block_seq, real_window_size, real_start, indexed_start,
                  ps[0:i].count('-'))

        if parent['start'] < 0:
            strand = -1
        else:
            strand = 1

        pid = percent_identity(block_seq, other['seq'][i:i + real_window_size])
        # Ignore 0% identity sequences
        if pid == 0:
            continue

        yield SeqFeature(
            FeatureLocation(
                real_start + indexed_start,
                real_start + indexed_start + nucl2prot * real_window_size
            ),
            type="match_part", strand=strand,
            qualifiers={
                "source": "progressiveMauve",
                'score': pid,
                # 'alignment': block_seq + '<br/>' + other['seq'][i:i + real_window_size],
                'qseq': block_seq,
                'sseq': other['seq'][i:i + real_window_size],
            }
        )


def reduce_subfeatures(subfeatures):
    prev_feature = None
    for feature in subfeatures:
        if prev_feature is None:
            prev_feature = feature
            continue

        if feature.location.start == prev_feature.location.end and prev_feature.qualifiers['score'] == feature.qualifiers['score']:
            prev_feature.location._end = feature.location.end
            # print prev_feature.location, feature.location, prev_feature.location.end, feature.location.start
            # pass
        else:
            yield prev_feature
            prev_feature = feature
    yield feature


def convert_xmfa_to_gff3(xmfa_file, sequences=None, window_size=1000, protein=False):
    label_convert = id_tn_dict(sequences)
    lcbs = parse_xmfa(xmfa_file)

    parent_records = {
        x: SeqRecord(Seq("ACTG", IUPAC.IUPACUnambiguousDNA), id=label_convert[x]['record_id'])
        for x in label_convert.keys()
    }

    for lcb_idx, lcb in enumerate(lcbs):
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
                    FeatureLocation(parent['start'], parent['end']),
                    type="match", strand=parent['strand'],
                    qualifiers={
                        'source': 'progressiveMauve',
                        'Target': label_convert[other['id']]['record_id'],
                        'Target_protein': other['comment'],
                        'ID': 'm_%s_%s_%s_%s' % (
                            lcb_idx, o_idx,
                            label_convert[parent['id']]['record_id'],
                            label_convert[other['id']]['record_id']
                        ),
                        'Name': other['comment'],
                    }
                )
                other_feature.sub_features = []

                subs = generate_subfeatures(parent, window_size, other, protein=protein)
                for subfeature in reduce_subfeatures(sorted(subs, key=lambda x: x.location.start)):
                    other_feature.sub_features.append(subfeature)

                parent_records[parent['id']].features.append(other_feature)

    for i in parent_records:
        yield parent_records[i]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert XMFA alignments to gff3', prog='xmfa2gff3')
    parser.add_argument('xmfa_file', type=argparse.FileType("r"), help='XMFA File')
    parser.add_argument('sequences', type=argparse.FileType("r"), help='Fasta file (in same order) passed to parent for reconstructing proper IDs')

    parser.add_argument('--window_size', type=int, help='Window size for analysis', default=1000)
    parser.add_argument('--protein', action='store_true', help='Protein XMFA file')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    try:
        os.makedirs('outdir')
    except:
        pass
    for result in convert_xmfa_to_gff3(**vars(args)):
        GFF.write([result], sys.stdout)
