#!/usr/bin/env python
import sys
import re
import hashlib
import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from collections import Counter
from phage_annotation_validator import excessive_gap, genes
import logging
logging.basicConfig(level=logging.DEBUG,  format='%(asctime)-15s %(message)s')
log = logging.getLogger('bemoan')


class GenomicUtils(object):
    """Genomic Analysis Utilities

    Refactored from the main bemoan software. All of this code is disgusting
    though, and really needs to be rewritten
    """

    @classmethod
    def normalize(cls, d, target=1.0):
        """Normalize all values in a dict
        """
        # http://stackoverflow.com/a/16418096
        raw = sum(d.values())
        factor = target/raw
        return {key: value * factor for key, value in d.iteritems()}

    @classmethod
    def find_orfs_with_trans2(cls, seq, trans_table=11, min_protein_length=20):
        seq_len = len(seq)
        # Leading M not needed because alternative starts
        pattern = re.compile('(^[A-Z]{3,})\*')
        start_codons = ['ATG', 'TTG', 'CTG', 'GTG', 'ATT', 'ATC', 'ATA']
        # For F and R
        for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
            # for each frame
            for frame in range(3):
                string_index = 0
                # Split into codons, groups of 3
                for codon in re.findall('..?.?', str(nuc[frame:])):
                    # Keep track of how far we've progressed through the
                    # reading frame
                    string_index+=3
                    # If it's a start codon
                    if codon in start_codons:
                        dna_seq = nuc[frame+string_index:]
                        # Trim to len % 3 == 0
                        if len(dna_seq) % 3 == 1:
                            dna_seq = dna_seq[0:-1]
                        elif len(dna_seq) % 3 == 2:
                            dna_seq = dna_seq[0:-2]
                        # Translate the entire frame
                        trans = str(dna_seq.translate(table=11))
                        for result in pattern.finditer(trans):  # Should there ever be more than one?
                            if len(result.group(1)) > min_protein_length:
                                if strand == 1:
                                    yield (
                                        string_index + frame,
                                        string_index + frame + 3 * len(result.group(1)),
                                        strand,
                                        result.group(1)
                                    )
                                else:
                                    yield (
                                        len(seq) - (string_index + frame + 3 * len(result.group(1))),
                                        len(seq) - (string_index + frame),
                                        strand,
                                        result.group(1)
                                    )

    @classmethod
    def features_with_bad_stops(cls, record, translation_table=11):
        bad_features = []
        for feature in record.features:
            protein = feature.extract(record).seq.translate(table=translation_table)
            if protein.count('*') > 1:
                bad_features.append(feature)
        return bad_features


class SDScore(object):
    """Class for naive SD scoring for a feature with upstream prepended
    """

    def __init__(self, sequence, upstream=300):
        """
            Initialize with a sequence with an added upstream portion, and the
            length of the upstream portion (used to recapitulate original
            start)
        """
        self.upstream = upstream
        self.sequence = sequence
        self.sds = [
            'aggaggt', 'ggaggt', 'aggagg', 'aggag', 'gaggt', 'ggagg', 'aggt',
            'gggt', 'gagg', 'gggg', 'agga', 'ggag', 'gga', 'gag', 'agg', 'ggt',
        ]

    def score_for_delta(self, delta, up_min=5, up_max=20):
        """
            Given a delta (new_start - original_start, negative is upstream),
            this will find the best SD in that sequence (if present, according
            to hardcoded list of SDs) and will return a "score"
        """
        seq_test = self.sequence[self.upstream + delta - up_max:self.upstream +
                                 delta - up_min]

        for i in self.sds:
            if i in seq_test.lower():
                return len(i)
        return 1


class PhageGeneCaller(object):
    """
        Manages application of gene annotations to a genome and generation of
        a genbank file
    """

    def __init__(self, genome=None):
        """Initialize with fasta file name
        """
        if genome is None:
            raise Exception("Must provide fasta genome file name")

        record = SeqIO.read(genome, 'fasta')
        self.record_id = record.id
        self.record = record
        self.record_len = len(record)

    def apply_annotations(self, gene_calls):
        """
            Given (usually) the results of CoalesceGeneCalls.coalesce, these
            gene calls are applied to the DNA sequence and returned

            This is an additive method
        """
        # Sorting to conform with bioperl behaviour
        for gene in sorted(gene_calls, key=lambda f: f.location.start):
            self.record.features.append(gene)

    def identify_empty_areas(self, allowed_overlap=30, min_open_frame=40):
        """Attempt to locate empty regions of the genome

        allowed_overlap is the distance which genes may overlap previously called genes

        min_open_frame is the smalled ORF to examine. Set to 50 AAs by default.
        """
        data = excessive_gap(self.record)
        return data[2]

    def annotate_empty_areas(self, region_list):
        """Attempt to locate empty regions of the genome, then do gene calling
        on them.

        Additionally this method blasts called features, and will colour them
        according to presence/absence of BLAST hits
        """
        for region in region_list:
            (start, end, unused) = region
            # Allow some overlap
            start -= 20
            end += 20
            # Fetch sequence
            sequence = self.record.seq[start:end]
            # Search in that sequence
            found_orfs = GenomicUtils.find_orfs_with_trans2(sequence, trans_table=11)
            log.info("Extracting region %s..%s. Found %s orfs", start, end, len(found_orfs))
            for (orf_start, orf_end, strand, protein) in found_orfs:
                # Region start, quit
                calc_start = start + orf_start
                calc_end = start + orf_end

                feature = SeqFeature(
                    FeatureLocation(calc_start, calc_end),
                    type="CDS",
                    strand=strand,
                    qualifiers={
                        'note': 'De-novo ORF called from six frame translation',
                        'color': '#aaffaa',
                    },
                )
                protein = str(feature.extract(self.record).seq.translate(table=11))

                if protein.count('*') > 1 and \
                        protein.index('*') < len(protein) - 1:
                    log.debug("Refusing to feature [%s %s %s] due to presence of %s stop codons" % (strand, calc_start, calc_end, protein.count('*')))
                else:
                    log.debug("Adding new feature [%s %s %s]" % (strand, calc_start, calc_end))

                print protein.count('*')

    def reannotate_existing_genes(self):
        orfs = []
        feature_map = {}
        for feature in self.record.features:
            try:
                protein = feature.extract(self.record).seq.translate(table=11)
                orfs.append([
                    feature.location.start,
                    feature.location.end,
                    0, 0,
                    feature.strand, protein, self.digest(str(protein))
                ])
                feature_map[self.digest(str(protein))] = feature
            except Exception, e:
                log.warn(e)
                log.info(feature)

        modifications = {}
        notes = {}

        for feature in self.record.features:
            protein = str(feature.extract(self.record).seq.translate(table=11))
            protein_hash = self.digest(protein)

            if protein_hash in modifications:
                # Here we create a temporary feature to check that if we were
                # to extend this protein according to our blast data, how the
                # results would look, would we start including some stop codons
                # by accident?
                tmp_feature = \
                    SeqFeature(FeatureLocation(feature.location.start.position
                                               + modifications[protein_hash],
                                               feature.location.end.position,
                                               strand=feature.location.strand),
                               type='domain')
                tmp_protein = str(tmp_feature.extract(self.record).seq.translate(table=11))
                if tmp_protein.count('*') > 1:
                    log.debug("Refusing to update start of %s from %s to %s due to presence of %s stop codons" %
                            (feature, feature.location.start.position, feature.location.start.position + modifications[protein_hash], protein.count('*')))
                else:
                    log.debug("Updating start of %s from %s to %s" % (feature, feature.location.start.position, feature.location.start.position + modifications[protein_hash]))

                    # Update before modifying feature location
                    msg = 'Updated start from %s to %s' % \
                        (feature.location.start.position,
                        feature.location.start.position +
                        modifications[protein_hash])
                    if 'note' in feature.qualifiers:
                        feature.qualifiers['note'].append(msg)
                    else:
                        feature.qualifiers['note'] = [msg]

                    feature.qualifiers['note'].append(str(notes[protein_hash]))

                    feature.location = FeatureLocation(
                        ExactPosition(feature.location.start.position +
                                    modifications[protein_hash]),
                        ExactPosition(feature.location.end.position),
                        feature.location.strand)

    def digest(self, seq):
        return hashlib.md5(seq).hexdigest()

    def _gen_seq_map(self, orf_list):
        seq_map = {}
        for (rs, re1, oq, oe, strand, protein) in orf_list:
            if crc not in seq_map:
                seq_map[crc] = [(rs, re1, oq, oe, strand)]
            else:
                seq_map[crc].append([rs, re1, oq, oe, strand])
        return seq_map

    def fix_bad_stops(self):
        for feature in GenomicUtils.features_with_bad_stops(self.record):
            # Find first stop
            first_stop = protein.find('*')
            strip_to_first_stop = len(protein) - first_stop

            log.debug("Updating stop of %s from %s to %s" % (feature,
                feature.location.start.position, feature.location.end.position
                - strip_to_first_stop))

            # Update before modifying feature location
            msg = 'Updated stop from %s to %s' % \
                (feature.location.end.position,
                 feature.location.end.position - strip_to_first_stop)

            self.append_qual(feature, 'note', msg)
            self.append_qual(feature, 'color', '4')

            feature.location = FeatureLocation(
                ExactPosition(feature.location.start.position),
                ExactPosition(feature.location.end.position - strip_to_first_stop),
                feature.location.strand)


    @classmethod
    def append_qual(cls, feature, qualifier, message):
        if 'qualifier' in feature.qualifiers:
            feature.qualifiers[qualifier].append(message)
        else:
            feature.qualifiers[qualifier] = [message]


class CoalesceGeneCalls(object):
    """
        Handles input of GFF data and coalescing of that data into a coherent
        set of annotations
    """

    def __init__(self, gff_data_sources=None):
        self.feature_groupings = {}
        self.features = []
        if gff_data_sources is not None:
            for source in gff_data_sources:
                self.add_gene_source(source)

    def add_gene_source(self, handle):
        """
            Add a set of features to the internal feature set. No
            calculations are run
        """
        for rec in GFF.parse(handle):
            for parent_feature in genes(rec.features, feature_type='CDS'):
                # Features are keyed on strand and end base, as that should be
                # shared amongst similar features
                if parent_feature.strand == 1:
                    fid = '%s:%s' % (parent_feature.location.end, parent_feature.strand)
                else:
                    fid = '%s:%s' % (parent_feature.location.start, parent_feature.strand)

                if fid not in self.feature_groupings:
                    self.feature_groupings[fid] = [parent_feature]
                else:
                    self.feature_groupings[fid].append(parent_feature)

    def coalesce(self):
        """
            Given the internal feature set (after loading of data), this
            function coalesces these into a single set of gene annotations
        """
        fixed_features = []
        for group_id in self.feature_groupings:
            log.debug('Processing %s' % group_id)
            feature_group = self.feature_groupings[group_id]
            fixed_features.append(
                self._process_group(feature_list=feature_group))
        log.info('Reduced to %s groups from %s features' %
                 (len(fixed_features), len(self.features)))
        return fixed_features

    @classmethod
    def _append_dict(cls, parent_dictionary, child_dictionary):
        for key in child_dictionary:
            if key not in parent_dictionary:
                parent_dictionary[key] = child_dictionary[key]
            else:
                parent_dictionary[key].extend(child_dictionary[key])
        return parent_dictionary

    @classmethod
    def _process_group(cls, feature_list=None):
        """
            Run the collapsing process.

            Currently this consists of just taking the longest feature, as that
            is vastly simpler. It doesn't make use of blast results which is
            what we'd really like. However, were we to blindly use blast
            results, we aren't gaining anything.

            To properly handle blast, we'd need to have a ranking of blast
            results, something to say "this phage is more relevant than that
            one"
        """
        end = feature_list[0].location.end
        strand = feature_list[0].location.strand
        start = None
        # We'll add another color if the start was updated
        changed_start = False
        qualifiers = {'note': []}
        # Pick Max Start/Longest ORF
        for feature in feature_list:
            if start is None:
                start = feature.location.start
            else:
                if feature.location.strand > 0:
                    if feature.location.start < start:
                        start = feature.location.start
                        changed_start = True
                else:
                    if feature.location.start > start:
                        start = feature.location.start
                        changed_start = True

            qualifiers = cls._append_dict(qualifiers, feature.qualifiers)
            for key in ['score', 'phase', 'ID']:
                if key in qualifiers:
                    del qualifiers[key]

        if changed_start:
            qualifiers['color'] = ['#ffaaaa']

        feature_loc = FeatureLocation(start, end, strand=strand)

        final_feature = SeqFeature(feature_loc, type="gene")
        final_cds = SeqFeature(feature_loc, type="CDS", qualifiers=qualifiers)
        final_feature.sub_features = [
            final_cds
        ]

        return final_feature


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='bemoan')
    parser.add_argument('genome', type=file, help='Genome Sequence')
    parser.add_argument('annotations', type=file, nargs='+', help='Parent GFF3 annotations')
    args = parser.parse_args()

    log.info("Loading Data")
    pgc = PhageGeneCaller(genome=args.genome)
    log.info("Gene Calls")
    gene_calls = CoalesceGeneCalls(gff_data_sources=args.annotations)
    log.info("Added annotations")
    pgc.apply_annotations(gene_calls.coalesce())
    # Annotate all empty streches first
    log.info("Annotation of empty ORFs")
    # We want to merge down the denovo features and the existing annotations.
    # Sometimes the de-novo caller will call multiple overlapping ORFs
    denovo_features = pgc.annotate_empty_areas(pgc.identify_empty_areas())
    #gene_calls2 = CoalesceGeneCalls()
    #gene_calls2.add_features(denovo_features)
    #pgc.apply_annotations(gene_calls2.coalesce())
    ## And then correct ALL the calls
    #log.info("Start Correction")
    #pgc.reannotate_existing_genes()
    #log.info("In-gene stop Corrections")
    #pgc.fix_bad_stops()

    GFF.write([pgc.record], sys.stdout)
