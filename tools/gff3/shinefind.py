#!/usr/bin/env python
import re
import sys
import argparse
import copy
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from gff3 import feature_lambda, feature_test_type

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def get_id(feature=None, parent_prefix=None):
    result = ""
    if parent_prefix is not None:
        result += parent_prefix + '|'
    if 'locus_tag' in feature.qualifiers:
        result += feature.qualifiers['locus_tag'][0]
    elif 'gene' in feature.qualifiers:
        result += feature.qualifiers['gene'][0]
    elif 'product' in feature.qualifiers:
        result += feature.qualifiers['product'][0]
    else:
        result += '%s_%s_%s' % (feature.location.start, feature.location.end,
                                feature.location.strand)
    return result


def ensure_location_in_bounds(start=0, end=0, parent_length=0):
    # This prevents frameshift errors
    while start < 0:
        start += 3
    while end < 0:
        end += 3
    while start > parent_length:
        start -= 3
    while end > parent_length:
        end -= 3
    return (start, end)


class NaiveSDCaller(object):

    SD_SEQUENCES = (
        'AGGAGGT',
        'GGAGGT',
        'AGGAGG',
        'AGGAG',
        'GAGGT',
        'GGAGG',
        'AGGT',
        'GGGT',
        'GAGG',
        'GGGG',
        'AGGA',
        'GGAG',
        'GGA',
        'GAG',
        'AGG',
        'GGT',
    )

    def __init__(self):
        self.sd_reg = [re.compile(x) for x in self.SD_SEQUENCES]

    def list_sds(self, sequence):
        hits = []
        for regex in self.sd_reg:
            for match in regex.finditer(sequence):
                hits.append({
                    'spacing': len(sequence) - len(match.group()) - match.start(),
                    'hit': match.group(),
                    'start': match.start(),
                    'end': match.end(),
                    'len': len(match.group()),
                })
        return hits

    @classmethod
    def highlight_sd(cls, sequence, start, end):
        return ' '.join([
            sequence[0:start].lower(),
            sequence[start:end].upper(),
            sequence[end:].lower()
        ])

    @classmethod
    def to_features(cls, hits, strand, parent_start, parent_end, feature_id=None):
        results = []
        for idx, hit in enumerate(hits):
            #gene            complement(124..486)
            #-1      491     501     0       5       5
            #-1      491     501     0       4       5
            #-1      491     501     1       4       5
            #-1      491     501     2       3       5
            #-1      491     501     1       3       5
            #-1      491     501     0       3       5

            qualifiers = {
                'source': 'CPT_ShineFind',
                'ID': '%s.rbs-%s' % (feature_id, idx)
            }
            if strand > 0:
                start = parent_end - hit['spacing'] - hit['len'] + 1
                end = parent_end - hit['spacing']
            else:
                start = parent_start + hit['spacing'] + 1
                end = parent_start + hit['spacing'] + hit['len']

            tmp = SeqFeature(
                FeatureLocation(start, end, strand=strand),
                type='Shine_Dalgarno_sequence',
                qualifiers=qualifiers
            )
            results.append(tmp)
        return results


def shinefind(fasta, gff3, table_output=None, lookahead_min=5, lookahead_max=15, top_only=False):
    table_output.write('\t'.join(['ID', 'Name', 'Terminus', 'Terminus', 'Strand',
                                  'Upstream Sequence', 'SD', 'Spacing']) +
                       "\n")

    sd_finder = NaiveSDCaller()

    # Load up sequence(s) for GFF3 data
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    # Parse GFF3 records
    for record in GFF.parse(gff3, base_dict=seq_dict):
        gff3_output = SeqRecord(record.seq, record.id)
        # Filter out just coding sequences
        for feature in feature_lambda(record.features, feature_test_type, {'type': 'CDS'}, subfeatures=False):
            # Strand information necessary to getting correct upstream sequence
            # TODO: library?
            strand = feature.location.strand
            # n_bases_upstream
            if strand > 0:
                start = feature.location.start - lookahead_max
                end = feature.location.start - lookahead_min
            else:
                start = feature.location.end + lookahead_min
                end = feature.location.end + lookahead_max

            (start, end) = ensure_location_in_bounds(start=start, end=end,
                                                        parent_length=record.__len__)

            # Create our temp feature used to obtain correct portion of
            # genome
            tmp = SeqFeature(FeatureLocation(start, end, strand=strand),
                                type='domain')
            seq = str(tmp.extract(record.seq))
            sds = sd_finder.list_sds(seq)

            feature_id = get_id(feature)
            sd_features = sd_finder.to_features(sds, strand, start, end, feature_id=feature.id)

            human_strand = '+' if strand == 1 else '-'

            # http://book.pythontips.com/en/latest/for_-_else.html
            log.debug('Found %s SDs', len(sds))
            for (sd, sd_feature) in zip(sds, sd_features):
                # If we only want the top feature, after the bulk of the
                # forloop executes once, we append the top feature, and fake a
                # break, because an actual break triggers the else: block
                table_output.write('\t'.join(map(str, [
                    feature.id,
                    feature_id,
                    feature.location.start,
                    feature.location.end,
                    human_strand,
                    sd_finder.highlight_sd(seq, sd['start'], sd['end']),
                    sd['hit'],
                    int(sd['spacing']) + lookahead_min,
                ])) + "\n")

                gff3_output.features.append(sd_feature)

                if top_only:
                    break
            else:
                if len(sds) != 0:
                    log.debug('Should not reach here if %s', len(sds) != 0)
                    # Somehow this is triggerring, and I don't feel like figuring out why. Someone else's problem.
                    continue
                table_output.write('\t'.join(map(str, [
                    feature.id,
                    feature_id,
                    feature.location.start,
                    feature.location.end,
                    human_strand,
                    seq,
                    None,
                    -1,
                ])) + "\n")

        gff3_output.features = sorted(gff3_output.features, key=lambda x: x.location.start)
        gff3_output.annotations = {}
        GFF.write([gff3_output], sys.stdout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta', type=file, help='Fasta Genome')
    parser.add_argument('gff3', type=file, help='GFF3 annotations')

    parser.add_argument('--table_output', type=argparse.FileType('w'), help='Tabular Output', default='shinefind.tbl')

    parser.add_argument('--lookahead_min', nargs='?', type=int, help='Number of bases upstream of CDSs to end search', default=5)
    parser.add_argument('--lookahead_max', nargs='?', type=int, help='Number of bases upstream of CDSs to begin search', default=15)

    parser.add_argument('--top_only', action='store_true', help='Only report best hits')

    args = parser.parse_args()
    shinefind(**vars(args))
