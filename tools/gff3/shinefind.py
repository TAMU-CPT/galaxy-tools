#!/usr/bin/env python
import re
import argparse
import copy
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def feature_lambda(feature_list, test, test_kwargs, subfeatures=True):
    # Either the top level set of [features] or the subfeature attribute
    for feature in feature_list:
        if test(feature, **test_kwargs):
            # Necessary? Or should we just wipe out actual feature's subfeatuers?
            if not subfeatures:
                feature_copy = copy.deepcopy(feature)
                feature_copy.sub_features = []
                yield feature_copy
            else:
                yield feature

        if hasattr(feature, 'sub_features'):
            for x in feature_lambda(feature.sub_features, test, test_kwargs, subfeatures=subfeatures):
                yield x


def feature_test(feature, **kwargs):
    return feature.type == kwargs['type']


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
    def to_features(cls, hits, strand, parent_start, parent_end):
        results = []
        for hit in hits:
            #gene            complement(124..486)
            #-1      491     501     0       5       5
            #-1      491     501     0       4       5
            #-1      491     501     1       4       5
            #-1      491     501     2       3       5
            #-1      491     501     1       3       5
            #-1      491     501     0       3       5

            if strand > 0:
                tmp = SeqFeature(FeatureLocation(
                            parent_end - hit['spacing'] - hit['len'] + 1,
                            parent_end - hit['spacing'],
                            strand=strand
                        ), type='RBS')
            else:
                tmp = SeqFeature(FeatureLocation(
                            parent_start + hit['spacing'] + 1,
                            parent_start + hit['spacing'] + hit['len'],
                            strand=strand
                        ), type='RBS')
            results.append(tmp)
        return results


def shinefind(fasta, gff3, gff3_output, table_output=None, lookahead_min=5, lookahead_max=15, top_only=False):
    table_output.write('\t'.join(['Name', 'Terminus', 'Terminus', 'Strand', 'Upstream Sequence', 'SD', 'Spacing']))

    sd_finder = NaiveSDCaller()

    # Load up sequence(s) for GFF3 data
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    # Parse GFF3 records
    for record in GFF.parse(gff3, base_dict=seq_dict):
        # Filter out just coding sequences
        for cds in feature_lambda(record.features, feature_test, {'type': 'CDS'}, subfeature=False):
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
            sd_features = sd_finder.to_features(sds, strand, start, end)

            feature_id = get_id(feature)
            human_strand = '+' if strand == 1 else '-'
            if len(sds) == 0:
                table_output.write('\t'.join([
                    feature_id,
                    feature.location.start,
                    feature.location.end,
                    human_strand,
                    seq,
                    None,
                    -1,
                ]))
                continue

            if top_only:
                table_output.write('\t'.join([
                    feature_id,
                    feature.location.start,
                    feature.location.end,
                    human_strand,
                    sd_finder.highlight_sd(seq, sds[0]['start'], sds[0]['end']),
                    sds[0]['hit'],
                    sds[0]['spacing'] + lookahead_min,
                ]))

                #gff3.append('\t'.join(map(str,[
                    #record.id,
                    #'CPT_ShineFind',
                    #'Shine_Dalgarno_sequence',
                    #sd_features[0].location.start,
                    #sd_features[0].location.end,
                    #'.',
                    #human_strand,
                    #'.',
                    #'ID=rbs_%s' % (feature_id, )
                #])))
            else:
                for sd in sds:
                    table_output.write('\t'.join([
                        feature_id,
                        feature.location.start,
                        feature.location.end,
                        human_strand,
                        sd_finder.highlight_sd(seq, sd['start'], sd['end']),
                        sd['hit'],
                        sd['spacing'] + lookahead_min,
                    ]))
                    #gff3.append('\t'.join(map(str,[
                        #record.id,
                        #'CPT_ShineFind',
                        #'Shine_Dalgarno_sequence',
                        #sd_features[0].location.start,
                        #sd_features[0].location.end,
                        #'.',
                        #human_strand,
                        #'.',
                        #'ID=rbs_%s' % (feature_id, )
                    #])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta', type=file, help='Fasta Genome')
    parser.add_argument('gff3', type=file, help='GFF3 annotations')

    parser.add_argument('gff3_output', help='GFF3 Output')
    parser.add_argument('--table_output', type=argparse.FileType('w'), help='Tabular Output', default='shinefind.tbl')

    parser.add_argument('--lookahead_min', nargs='?', type=int, help='Number of bases upstream of CDSs to end search', default=5)
    parser.add_argument('--lookahead_max', nargs='?', type=int, help='Number of bases upstream of CDSs to begin search', default=15)

    parser.add_argument('--top_only', action='store_true', help='Only report best hits')

    args = parser.parse_args()
    shinefind(**vars(args))
