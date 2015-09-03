#!/usr/bin/env python
import sys
import re
import hashlib
import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from collections import Counter

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
        answer = []
        seq_len = len(seq)
        # Leading M not needed because alternative starts
        pattern = re.compile('(^[A-Z]{3,})\*')
        start_codons = ['ATG', 'TTG', 'CTG', 'GTG']
        # Table 11 says these are starts too? 'ATT', 'ATC', 'ATA'
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
                                    answer.append((
                                        # Index to get our left shift, in addition to frame.
                                        string_index + frame - 3,
                                        # Then we take 3*end (length
                                        # basically), and add 3 to account for
                                        # added stop codon
                                        3*result.end(1) + 3,
                                        strand,
                                        result.group(1)
                                    ))
                                else:
                                    log.debug(result.group(1))
                                    log.debug(' '.join(map(str, [result.start(1),
                                                                 result.end(1),
                                                                 3*result.start(1),
                                                                 3*result.end(1),
                                                                 frame,
                                                                 string_index])))
                                    answer.append((
                                        # See below, then subtract length of protein found on top of that
                                        seq_len - string_index - 3*result.end(1) - frame + 3,
                                        # Here we have the "end" (right most
                                        # portion), which takes length of input
                                        # sequence, and subtracts our progress
                                        # along that (since we're looking at
                                        # that sequence backwards)
                                        3*result.end(1) + 3,
                                        strand,
                                        result.group(1)
                                    ))
        answer.sort()
        return answer

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

    def get_annotations(self):
        """Returns annotations currently applied to genome object"""
        return self.record.features

    def clear_annotations(self):
        """Removes all annotations from object
        """
        self.record.features = []

    def identify_empty_areas(self, allowed_overlap=30, min_open_frame=40):
        """Attempt to locate empty regions of the genome

        allowed_overlap is the distance which genes may overlap previously called genes

        min_open_frame is the smalled ORF to examine. Set to 50 AAs by default.
        """
        annotated = set(range(len(self.record)))
        for feature in self.record.features:
            annotated = annotated - set(range(feature.location.start,
                                              feature.location.end))
        region_list = []
        current_region_start = None
        last_nt = None

        for nt in sorted(list(annotated)):
            if current_region_start is None:
                current_region_start = nt
            if last_nt is None:
                last_nt = nt
            if nt != last_nt + 1:
                region_list.append([current_region_start - allowed_overlap,
                                    last_nt + allowed_overlap])
                last_nt = nt
                current_region_start = nt
            else:
                last_nt = nt
        region_list.append([current_region_start - allowed_overlap, last_nt +
                            allowed_overlap])

        # remove short empty regions
        region_list = [l for l in region_list if abs(l[1]-l[0]) >
                       min_open_frame*3]

        return region_list

    def annotate_empty_areas(self, region_list):
        """Attempt to locate empty regions of the genome, then do gene calling
        on them.

        Additionally this method blasts called features, and will colour them
        according to presence/absence of BLAST hits
        """
        new_possible_orfs = []
        for region in region_list:
            (start, end) = region
            log.info("\n\nExtracting region %s..%s" % (start, end))
            sequence = self.record.seq[start:end]
            orf_list = GenomicUtils.find_orfs_with_trans2(sequence, trans_table=11)
            for (orf_start, orf_end, strand, protein) in orf_list:
                new_possible_orfs.append(( start, end, orf_start,
                                          orf_end, strand, str(protein),
                                          self.crc32(str(protein)) ))

        blast_stdout = self._blast_orfs(new_possible_orfs)
        seq_map = self._gen_seq_map(new_possible_orfs)

        # Find all listed protein hashes in the blast results, we're only
        # interested in de-novo features if they have blast results. Raw orf
        # calling isn't likely to be succesfully without evidence
        good_ids = {}
        for line in blast_stdout.split('\n'):
            if len(line) > 1:
                good_ids[line.split('\t')[0]] = True
        log.debug("Found these blast hits: " + str(good_ids))

        new_blast_features = []
        # For all of the protein hashes found in blast results, if they're in
        # our sequence map (and they bloody well should be, but it happens?
        # Sometimes?)
        for crc_id in [x for x in seq_map]:
            for feature_info in seq_map[crc_id]:
                # Region start, quit
                (rs, rq, dist_from_left_end, seq_length, strand) = feature_info
                calc_start = rs + dist_from_left_end
                calc_end = rs + dist_from_left_end + seq_length

                feature = SeqFeature(FeatureLocation(calc_start, calc_end),
                                     type="CDS", strand=strand)
                protein = str(feature.extract(self.record).seq.translate(table=11))

                if crc_id in good_ids:
                    # Has blast data backing it up
                    self.append_qual(feature, 'note', 'De-novo ORF called with BLAST hit')
                    self.append_qual(feature, 'color', '3')
                    # Translate to check for stops
                    # This actually catches some
                    if protein.count('*') > 1 and \
                            protein.index('*') < len(protein) - 1:
                        log.debug("Refusing to add BLAST de-novo protein [%s %s %s] due to presence of %s stop codons" % (strand, calc_start, calc_end, protein.count('*')))
                    else:
                        log.debug("Adding BLAST de-novo protein [%s %s %s]" % (strand, calc_start, calc_end))
                        new_blast_features.append(feature)
                else:
                    self.append_qual(feature, 'note', 'De-novo ORF called without BLAST hit')
                    self.append_qual(feature, 'color', '2')
                    # DON'T APPEND FEATURE


        return new_blast_features

    def analyse_group(self, blast_hit_list, sequence_info):
        # Ignore empty lists, no blast hits shouldn't get here but just in case
        if len(blast_hit_list) == 0:
            return None

        # Check that we're completely aligned:
        alignment = {
            'ends_exact': 0,     # Ends exactly matched, we assume that we aren't
                                 # missing stop codons/gene callers aren't missing
                                 # stop codons either
            'start_exact': 0,    # Starts exactly matched
            'longer_than_db': 0,   # query is longer than hit, maybe move downstream
            'shorter_than_db': 0,  # query is shorter, maybe move upstream
            'diffs': [],         # differences in starts
        }
        for hit in blast_hit_list:
            if hit.qend == hit.send:
                alignment['ends_exact'] += 1
            if hit.qstart == hit.sstart:
                alignment['start_exact'] += 1
            elif hit.qstart < hit.sstart:
                alignment['shorter_than_db'] += 1
            elif hit.qstart > hit.sstart:
                alignment['longer_than_db'] += 1
            alignment['diffs'].append(int(hit.qstart) - int(hit.sstart))

        # We've got a set of alignments. We're under the unique constraint of
        # wanting longer alignments if they exist, thus we bias towards longest
        # possible gene length. We'll aim for a reasonable % of the hits being
        # longer than our query sequence before we change the start.

        # This method no longer appears to bias towards upstream hits

        # If the number of hits in the DB that are longer than our query
        # account for >15% of the hits, then we're interested.
        hit_hist = {k: v for k, v in
                    GenomicUtils.normalize(Counter(alignment['diffs'])).iteritems() if
                    v > .10}
        if len(hit_hist) > 0:
            # We have some possibilities. We might have multiple ones.
            # We have some k/v pairs in hit_hist representing a start change
            # with at least 10% of hits
            sd_score = SDScore(sequence_info['seq'],
                               upstream=sequence_info['upstream_distance'])
            for key in hit_hist:
                # We'll need to evaluate the SD sequence available to that
                # start
                best = sd_score.score_for_delta(key)
                # Multiply this with our fractional hit number...
                hit_hist[key] *= best

            # Then return the highest scoring one!
            #inverse = [(value, key) for key, value in hit_hist.items()]
            #best_upstream_start = max(inverse)[1]
            #return best_upstream_start
            return hit_hist
        else:
            return None

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
                    feature.strand, protein, self.crc32(str(protein))
                ])
                feature_map[self.crc32(str(protein))] = feature
            except Exception, e:
                log.warn(e)
                log.info(feature)

        modifications = {}
        notes = {}

        for feature in self.record.features:
            protein = str(feature.extract(self.record).seq.translate(table=11))
            protein_hash = self.crc32(protein)

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

    def crc32(self, seq):
        return hashlib.md5(seq).hexdigest()
        #return '%08x' % (binascii.crc32(seq) & 0xffffffff)

    def _gen_seq_map(self, orf_list):
        seq_map = {}
        for (rs, re1, oq, oe, strand, protein, crc) in orf_list:
            if crc not in seq_map:
                seq_map[crc] = [[rs, re1, oq, oe, strand]]
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

    def append_qual(self, feature, qualifier, message):
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
            # TODO: Provide method for obtaining RBS features/other types
            for feature in [f for f in rec.features if f.type == 'CDS']:
                self.features.append(feature)

    def add_features(self, feature_list):
        for feature in feature_list:
            if feature.type == "CDS":
                self.features.append(feature)

    def coalesce(self):
        """
            Given the internal feature set (after loading of data), this
            function coalesces these into a single set of gene annotations
        """

        for feature in self.features:
            # We assume that all features have correct stop data
            if feature.strand == 1:
                fid = '%s:%s' % (feature.location.end, feature.strand)
            else:
                fid = '%s:%s' % (feature.location.start, feature.strand)
            if fid not in self.feature_groupings:
                self.feature_groupings[fid] = [feature]
            else:
                self.feature_groupings[fid].append(feature)

        if len(self.features) == 34:
            import pprint; pprint.pprint(self.feature_groupings)

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
        if not isinstance(feature_list, list) and feature_list is not None:
            raise Exception("Must provide list of features")
        end = feature_list[0].location.end
        strand = feature_list[0].location.strand
        start = None
        # We'll add another color if the start was updated
        changed_start = False
        qualifiers = {'note':[]}
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

            for key in feature.qualifiers:
                for value in feature.qualifiers[key]:
                    if key not in qualifiers:
                        qualifiers[key] = []

                    if key in ['score', 'phase', 'ID']:
                        # Ignore completely
                        pass
                    elif key is 'source':
                        qualifiers['note'].append('%s..%s Source %s' % (start, end, value))
                    elif key not in ['color', 'colour']:
                        qualifiers[key].append('%s..%s %s' % (start, end, value))
                    else:
                        qualifiers[key] = value
        if changed_start:
            qualifiers['color'] = ['4']
        final_feature = SeqFeature(FeatureLocation(start, end, strand=strand),
                                   type="CDS", qualifiers=qualifiers)
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
    empty_region_list = pgc.identify_empty_areas()
    denovo_features = pgc.annotate_empty_areas(empty_region_list)
    gene_calls2 = CoalesceGeneCalls()
    gene_calls2.add_features(denovo_features)
    pgc.apply_annotations(gene_calls2.coalesce())
    # And then correct ALL the calls
    log.info("Start Correction")
    pgc.reannotate_existing_genes()
    log.info("In-gene stop Corrections")
    pgc.fix_bad_stops()

    GFF.write([pgc.record], sys.stdout)
