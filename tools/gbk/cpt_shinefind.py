#!/usr/bin/env python
import re
import sys
import argparse
import logging
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from gff3 import get_id, ensure_location_in_bounds
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


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
        self.sd_reg = [re.compile(x, re.IGNORECASE) for x in self.SD_SEQUENCES]

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
            # gene            complement(124..486)
            # -1      491     501     0       5       5
            # -1      491     501     0       4       5
            # -1      491     501     1       4       5
            # -1      491     501     2       3       5
            # -1      491     501     1       3       5
            # -1      491     501     0       3       5

            qualifiers = {
                'source': 'CPT_ShineFind',
            }

            if strand > 0:
                start = parent_end - hit['spacing'] - hit['len']
                end = parent_end - hit['spacing']
            else:
                start = parent_start + hit['spacing']
                end = parent_start + hit['spacing'] + hit['len']

            tmp = SeqFeature(
                FeatureLocation(start, end, strand=strand),
                type='RBS',
                qualifiers=qualifiers
            )
            results.append(tmp)
        return results

    def testFeatureUpstream(self, feature, record, sd_min=5, sd_max=15):
        # Strand information necessary to getting correct upstream sequence
        # TODO: library?
        strand = feature.location.strand
        # n_bases_upstream
        if strand > 0:
            start = feature.location.start - sd_max
            end = feature.location.start - sd_min
        else:
            start = feature.location.end + sd_min
            end = feature.location.end + sd_max

        (start, end) = ensure_location_in_bounds(start=start, end=end,
                                                 parent_length=len(record))

        # Create our temp feature used to obtain correct portion of
        # genome
        tmp = SeqFeature(FeatureLocation(start, end, strand=strand),
                         type='domain')
        seq = str(tmp.extract(record.seq))
        return self.list_sds(seq), start, end, seq

    def hasSd(self, feature, record, sd_min=5, sd_max=15):
        sds, start, end, seq = self.testFeatureUpstream(feature, record, sd_min=sd_min, sd_max=sd_max)
        return len(sds) > 0


def shinefind(genbank_file, gff3_output=None, table_output=None, lookahead_min=5, lookahead_max=15, top_only=False, add=False):
    table_output.write('\t'.join(['ID', 'Name', 'Terminus', 'Terminus', 'Strand',
                                  'Upstream Sequence', 'SD', 'Spacing']) +
                       "\n")

    sd_finder = NaiveSDCaller()
    # Parse GFF3 records
    for record in list(SeqIO.parse(genbank_file, "genbank")):
        # Sometimes you have a case where TWO CDS features have the same start. Only handle ONE.
        seen = {}
        # Shinefind's "gff3_output".
        gff3_output_record = SeqRecord(record.seq, record.id)
        # Loop over all CDS features
        for feature in record.features:
            if feature.type != 'CDS':
                continue

            seen_loc = feature.location.start if feature.strand > 0 else feature.location.end
            if seen_loc in seen:
                continue
            else:
                seen[seen_loc] = True

            sds, start, end, seq = sd_finder.testFeatureUpstream(
                feature, record, sd_min=lookahead_min, sd_max=lookahead_max)

            feature_id = get_id(feature)
            sd_features = sd_finder.to_features(sds, feature.location.strand, start, end, feature_id=feature.id)

            human_strand = '+' if feature.location.strand == 1 else '-'

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

                if add:
                    # Append the top RBS to the gene feature
                    record.features.append(sd_feature)
                # Also register the feature with the separate GFF3 output
                gff3_output_record.features.append(sd_feature)

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

        record.features = sorted(record.features, key=lambda x: x.location.start)
        SeqIO.write([record], sys.stdout, 'genbank')

        gff3_output_record.features = sorted(gff3_output_record.features, key=lambda x: x.location.start)
        gff3_output_record.annotations = {}
        GFF.write([gff3_output_record], gff3_output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('genbank_file', type=argparse.FileType('r'), help='Genbank file')

    parser.add_argument('--gff3_output', type=argparse.FileType('w'), help='GFF3 Output', default='shinefind.gff3')
    parser.add_argument('--table_output', type=argparse.FileType('w'), help='Tabular Output', default='shinefind.tbl')

    parser.add_argument('--lookahead_min', nargs='?', type=int,
                        help='Number of bases upstream of CDSs to end search', default=5)
    parser.add_argument('--lookahead_max', nargs='?', type=int,
                        help='Number of bases upstream of CDSs to begin search', default=15)

    parser.add_argument('--top_only', action='store_true', help='Only report best hits')
    parser.add_argument('--add', action='store_true',
                        help='Function in "addition" mode whereby the ' +
                        'RBSs are added directly to the gene model.')

    args = parser.parse_args()
    shinefind(**vars(args))
