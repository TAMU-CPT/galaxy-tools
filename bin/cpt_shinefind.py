#!/usr/bin/env python
import sys
import re
import logging
logging.basicConfig(level=logging.INFO)
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import StringIO

def get_id(feature=None, parent_prefix=None, idx=None):
    result = ""
    if parent_prefix is not None:
        result += parent_prefix + '|'
    if idx is not None:
        result += '%03d_' % idx
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
                })
        return hits

    @classmethod
    def highlight_sd(cls, sequence, start, end):
        return ' '.join([
            sequence[0:start].lower(),
            sequence[start:end].upper(),
            sequence[end:].lower()
        ])

def shinefind(genbank_file, table_output, gff3_output, lookahead_min=5, lookahead_max=15, top_only=False):
    output = StringIO.StringIO()
    records = list(SeqIO.parse(genbank_file, "genbank"))

    results = [['Name', 'Terminus', 'Terminus', 'Strand', 'Upstream Sequence', 'SD', 'Spacing']]
    gff3 = ['##gff-version 3']

    sd_finder = NaiveSDCaller()

    for record in records:
        for feature in record.features:
            if feature.type in ('CDS',):
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
                human_strand = '+' if strand == 1 else '-'
                if len(sds) == 0:
                    results.append([
                        feature_id,
                        feature.location.start,
                        feature.location.end,
                        human_strand,
                        seq,
                        None,
                        -1,
                    ])
                    continue

                if top_only:
                    results.append([
                        feature_id,
                        feature.location.start,
                        feature.location.end,
                        human_strand,
                        sd_finder.highlight_sd(seq, sds[0]['start'], sds[0]['end']),
                        sds[0]['hit'],
                        sds[0]['spacing'] + lookahead_min,
                    ])

                    gff3.append('\t'.join(map(str,[
                        record.id,
                        'CPT_ShineFind',
                        'Shine_Dalgarno_sequence',
                        sds[0]['hit'],
                        sds[0]['spacing'] + lookahead_min,
                        '.',
                        human_strand,
                        '.',
                        'ID=rbs_%s' % (feature_id, )
                    ])))
                else:
                    for sd in sds:
                        results.append([
                            feature_id,
                            feature.location.start,
                            feature.location.end,
                            human_strand,
                            sd_finder.highlight_sd(seq, sd['start'], sd['end']),
                            sd['hit'],
                            sd['spacing'] + lookahead_min,
                        ])
                        gff3.append('\t'.join(map(str,[
                            record.id,
                            'CPT_ShineFind',
                            'Shine_Dalgarno_sequence',
                            sd['hit'],
                            sd['spacing'] + lookahead_min,
                            '.',
                            human_strand,
                            '.',
                            'ID=rbs_%s' % (feature_id, )
                        ])))

    human_table = '\n'.join(['\t'.join(map(str, row)) for row in results])
    with open(table_output, 'w') as handle:
        handle.write(human_table)

    with open(gff3_output, 'w') as handle:
        handle.write('\n'.join(gff3))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('genbank_file', type=file, nargs='?',
                                        help='Genbank file')
    parser.add_argument('table_output', help='Tabular Output')
    parser.add_argument('gff3_output', help='GFF3 Output')
    parser.add_argument('--lookahead_min', nargs='?', type=int, help='Number of bases upstream of CDSs to end search', default=5)
    parser.add_argument('--lookahead_max', nargs='?', type=int, help='Number of bases upstream of CDSs to begin search', default=15)

    parser.add_argument('--top_only', action='store_true', help='Only report best hits')

    args = parser.parse_args()
    shinefind(**vars(args))
