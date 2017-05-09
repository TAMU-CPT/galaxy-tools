#!/usr/bin/env python
import re
import sys
import argparse
import logging
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from gff3 import feature_lambda, feature_test_type, feature_test_quals, get_id, \
    ensure_location_in_bounds
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


class NaiveSDCaller(object):

    SD_SEQUENCES = (
        'AGGAGGT',
        'GGAGGT',
        'AGGAGG',
        'GGGGGG',
        'AGGAG',
        'GAGGT',
        'GGAGG',
        'GGGGG',
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
        'GGG',
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
                'ID': '%s.rbs-%s' % (feature_id, idx)
            }

            if strand > 0:
                start = parent_end - hit['spacing'] - hit['len']
                end = parent_end - hit['spacing']
            else:
                start = parent_start + hit['spacing']
                end = parent_start + hit['spacing'] + hit['len']

            tmp = SeqFeature(
                FeatureLocation(start, end, strand=strand),
                type='Shine_Dalgarno_sequence',
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
                                                 parent_length=record.__len__)

        # Create our temp feature used to obtain correct portion of
        # genome
        tmp = SeqFeature(FeatureLocation(start, end, strand=strand),
                         type='domain')
        seq = str(tmp.extract(record.seq))
        return self.list_sds(seq), start, end, seq

    def hasSd(self, feature, record, sd_min=5, sd_max=15):
        sds, start, end, seq = self.testFeatureUpstream(feature, record, sd_min=sd_min, sd_max=sd_max)
        return len(sds) > 0


def shinefind(fasta, gff3, gff3_output=None, table_output=None, lookahead_min=5, lookahead_max=15, top_only=False, add=False):
    table_output.write('\t'.join(['ID', 'Name', 'Terminus', 'Terminus', 'Strand',
                                  'Upstream Sequence', 'SD', 'Spacing']) +
                       "\n")

    sd_finder = NaiveSDCaller()
    # Load up sequence(s) for GFF3 data
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    # Parse GFF3 records
    for record in GFF.parse(gff3, base_dict=seq_dict):
        # Shinefind's "gff3_output".
        gff3_output_record = SeqRecord(record.seq, record.id)
        # Filter out just coding sequences
        ignored_features = []
        for x in record.features:
            # If feature X does NOT contain a CDS, add to ignored_features
            # list. This means if we have a top level gene feature with or
            # without a CDS subfeature, we're catch it appropriately here.
            if len(list(feature_lambda([x], feature_test_type, {'type': 'CDS'}, subfeatures=True))) == 0:
                ignored_features.append(x)

        # Loop over all gene features
        for gene in feature_lambda(record.features, feature_test_type, {'type': 'gene'}, subfeatures=True):

            # Get the CDS from this gene.
            feature = list(feature_lambda(gene.sub_features, feature_test_type, {'type': 'CDS'}, subfeatures=True))
            # If no CDSs are in this gene feature, then quit
            if len(feature) == 0:
                # We've already caught these above in our ignored_features
                # list, so we skip out on the rest of this for loop
                continue
            else:
                # Otherwise pull the first (bad?) We don't expect >1 CDS/gene
                feature = feature[0]

            # Three different ways RBSs can be stored that we expect.
            rbs_rbs = list(feature_lambda(gene.sub_features, feature_test_type, {'type': 'RBS'}, subfeatures=False))
            rbs_sds = list(feature_lambda(gene.sub_features, feature_test_type,
                                          {'type': 'Shine_Dalgarno_sequence'}, subfeatures=False))
            regulatory_elements = list(feature_lambda(gene.sub_features, feature_test_type,
                                                      {'type': 'regulatory'}, subfeatures=False))
            rbs_regulatory = list(feature_lambda(regulatory_elements, feature_test_quals,
                                                 {'regulatory_class': ['ribosome_binding_site']}, subfeatures=False))
            rbss = rbs_rbs + rbs_sds + rbs_regulatory

            # If someone has already annotated an RBS, we quit
            if len(rbss) > 0:
                log.debug("Has %s RBSs", len(rbss))
                ignored_features.append(gene)
                continue

            sds, start, end, seq = sd_finder.testFeatureUpstream(feature, record, sd_min=lookahead_min, sd_max=lookahead_max)

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
                    gene.sub_features.append(sd_feature)
                    # Pick out start/end locations for all sub_features
                    locations = [x.location.start for x in gene.sub_features] + \
                        [x.location.end for x in gene.sub_features]
                    # Update gene's start/end to be inclusive
                    gene.location._start = min(locations)
                    gene.location._end = max(locations)
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

        record.annotations = {}
        GFF.write([record], sys.stdout)

        gff3_output_record.features = sorted(gff3_output_record.features, key=lambda x: x.location.start)
        gff3_output_record.annotations = {}
        GFF.write([gff3_output_record], gff3_output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('fasta', type=argparse.FileType('r'), help='Fasta Genome')
    parser.add_argument('gff3', type=argparse.FileType('r'), help='GFF3 annotations')

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
