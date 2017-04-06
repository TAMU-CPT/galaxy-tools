#!/usr/bin/env python
# vim: set fileencoding=utf-8
import os
import json
import math
import argparse
import itertools
import logging
from gff3 import feature_lambda, \
    coding_genes, genes, get_gff3_id, feature_test_location, get_rbs_from, nice_name
from shinefind import NaiveSDCaller
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from jinja2 import Environment, FileSystemLoader
from cpt import MGAFinder
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name='pav')

# Path to script, required because of Galaxy.
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
# Path to the HTML template for the report

ENCOURAGEMENT = (
    (100, 'Perfection itself!'),
    (90, 'Amazing!'),
    (80, 'Not too bad, a few minor things to fix...'),
    (70, 'Some issues to address'),
    (50, """Issues detected! </p><p class="text-muted">Have you heard of the
     <a href="https://cpt.tamu.edu">CPT</a>\'s Automated Phage Annotation
     Pipeline?"""),
    (0, """<b>MAJOR</b> issues detected! Please consider using the
     <a href="https://cpt.tamu.edu">CPT</a>\'s Automated Phage Annotation Pipeline"""),
)


def gen_qc_feature(start, end, message, strand=0, id_src=None):
    kwargs = {
        'qualifiers': {
            'note': [message]
        }
    }
    if id_src is not None:
        kwargs['id'] = id_src.id
        kwargs['qualifiers']['Name'] = id_src.qualifiers.get('Name', [])

    return SeqFeature(
        FeatureLocation(start, end, strand=strand),
        **kwargs
    )


def __ensure_location_in_bounds(start=0, end=0, parent_length=0):
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


def missing_rbs(record, lookahead_min=5, lookahead_max=15):
    """
    Identify gene features with missing RBSs

    This "looks ahead" 5-15 bases ahead of each gene feature, and checks if
    there's an RBS feature in those bounds.

    The returned data is a set of genes with the RBS sequence in the __upstream
    attribute, and a message in the __message attribute.
    """
    results = []
    good = 0
    bad = 0
    qc_features = []
    sd_finder = NaiveSDCaller()

    any_rbss = False

    for gene in coding_genes(record.features):
        # Check if there are RBSs, TODO: make this recursive. Each feature in
        # gene.sub_features can also have sub_features.
        rbss = get_rbs_from(gene)
        # No RBS found
        if len(rbss) == 0:
            # Get the sequence lookahead_min to lookahead_max upstream
            if gene.strand > 0:
                start = gene.location.start - lookahead_max
                end = gene.location.start - lookahead_min
            else:
                start = gene.location.end + lookahead_min
                end = gene.location.end + lookahead_max
            # We have to ensure the feature is ON the genome, otherwise we may
            # be trying to access a location outside of the length of the
            # genome, which would be bad.
            (start, end) = __ensure_location_in_bounds(start=start, end=end,
                                                       parent_length=record.__len__)
            # Temporary feature to extract sequence
            tmp = SeqFeature(FeatureLocation(start, end, strand=gene.strand),
                             type='domain')
            # Get the sequence
            seq = str(tmp.extract(record.seq))
            # Set the default properties
            gene.__upstream = seq.lower()
            gene.__message = "No RBS annotated, None found"

            # Try and do an automated shinefind call
            sds = sd_finder.list_sds(seq)
            if len(sds) > 0:
                sd = sds[0]
                gene.__upstream = sd_finder.highlight_sd(seq.lower(), sd['start'], sd['end'])
                gene.__message = "Unannotated but valid RBS"

            qc_features.append(gen_qc_feature(start, end, 'Missing RBS', strand=gene.strand, id_src=gene))

            bad += 1
            results.append(gene)
        else:
            if len(rbss) > 1:
                log.warn("%s RBSs found for gene %s", rbss[0].id, get_gff3_id(gene))
            any_rbss = True
            # get first RBS/CDS
            cds = list(genes(gene.sub_features, feature_type='CDS'))[0]
            rbs = rbss[0]

            # Get the distance between the two
            if gene.strand > 0:
                distance = cds.location.start - rbs.location.end
            else:
                distance = rbs.location.start - cds.location.end

            # If the RBS is too far away, annotate that
            if distance > lookahead_max:
                gene.__message = "RBS too far away (%s nt)" % distance

                qc_features.append(gen_qc_feature(
                    rbs.location.start,
                    rbs.location.end,
                    gene.__message,
                    strand=gene.strand,
                    id_src=gene,
                ))

                bad += 1
                results.append(gene)
            else:
                good += 1

    return good, bad, results, qc_features, any_rbss

# modified from get_orfs_or_cdss.py
# -----------------------------------------------------------


def require_sd(data, record, chrom_start, sd_min, sd_max):
    sd_finder = NaiveSDCaller()
    for putative_gene in data:
        if putative_gene[2] > 0:  # strand
            start = chrom_start + putative_gene[0] - sd_max
            end = chrom_start + putative_gene[0] - sd_min
        else:
            start = chrom_start + putative_gene[1] + sd_min
            end = chrom_start + putative_gene[1] + sd_max

        (start, end) = __ensure_location_in_bounds(start=start, end=end,
                                                   parent_length=record.__len__)
        tmp = SeqFeature(FeatureLocation(
            start, end, strand=putative_gene[2]), type='domain')
        # Get the sequence
        seq = str(tmp.extract(record.seq))
        sds = sd_finder.list_sds(seq)
        if len(sds) > 0:
            yield putative_gene + (start, end)


def excessive_gap(record, excess=50, excess_divergent=200, min_gene=30, slop=30, lookahead_min=5, lookahead_max=15):
    """
    Identify excessive gaps between gene features.

    Default "excessive" gap size is 10, but that should likely be larger.
    """
    results = []
    good = 0
    bad = 0

    contiguous_regions = []

    sorted_genes = sorted(genes(record.features), key=lambda feature: feature.location.start)
    if len(sorted_genes) == 0:
        log.warn("NO GENES FOUND")
        return good, bad, results, []

    current_gene = None
    for gene in sorted_genes:
        # If the gene's start is contiguous to the "current_gene", then we
        # extend current_gene
        log.debug('gene.id', gene.id)
        for cds in genes(gene.sub_features, feature_type='CDS'):
            log.debug('\t%s %s', cds.id, cds.location)
            if current_gene is None:
                current_gene = [
                    int(cds.location.start),
                    int(cds.location.end)
                ]

            if cds.location.start <= current_gene[1] + excess:
                # Don't want to decrease size
                if int(cds.location.end) >= current_gene[1]:
                    current_gene[1] = int(cds.location.end)
            else:
                # If it's discontiguous, we append the region and clear.
                contiguous_regions.append(current_gene)
                current_gene = [int(cds.location.start), int(cds.location.end)]

    # This generally expected that annotations would NOT continue unto the end
    # of the genome, however that's a bug, and we can make it here with an
    # empty contiguous_regions list
    contiguous_regions.append(current_gene)

    for i in range(len(contiguous_regions) + 1):
        if i == 0:
            a = (1, 1)
            b = contiguous_regions[i]
        elif i >= len(contiguous_regions):
            a = contiguous_regions[i - 1]
            b = (len(record.seq), None)
        else:
            a = contiguous_regions[i - 1]
            b = contiguous_regions[i]

        gap_size = abs(b[0] - a[1])
        if gap_size > min(excess, excess_divergent):
            a_feat_l = itertools.islice(feature_lambda(sorted_genes, feature_test_location, {'loc': a[1]}, subfeatures=False), 1)
            b_feat_l = itertools.islice(feature_lambda(sorted_genes, feature_test_location, {'loc': b[0]}, subfeatures=False), 1)

            try:
                a_feat = next(a_feat_l)
            except StopIteration:
                # Triggers on end of genome
                a_feat = None
            try:
                b_feat = next(b_feat_l)
            except StopIteration:
                # Triggers on end of genome
                b_feat = None

            result_obj = [
                a[1],
                b[0],
                None if not a_feat else a_feat.location.strand,
                None if not b_feat else b_feat.location.strand
            ]

            if a_feat is None or b_feat is None:
                if gap_size > excess_divergent:
                    results.append(result_obj)
            else:
                if a_feat.location.strand == b_feat.location.strand and gap_size > excess:
                    results.append(result_obj)
                elif a_feat.location.strand != b_feat.location.strand and gap_size > excess_divergent:
                    results.append(result_obj)

    better_results = []
    qc_features = []
    of = MGAFinder(11, 'CDS', 'closed', min_gene)
    # of = OrfFinder(11, 'CDS', 'closed', min_gene)

    for result_obj in results:
        start = result_obj[0]
        end = result_obj[1]
        f = gen_qc_feature(start, end, 'Excessive gap, %s bases' % abs(end - start))
        qc_features.append(f)
        putative_genes = of.putative_genes_in_sequence(str(record[start - slop:end + slop].seq))
        putative_genes = list(require_sd(putative_genes, record, start, lookahead_min, lookahead_max))
        for putative_gene in putative_genes:
            # (0, 33, 1, 'ATTATTTTATCAAAACGCTTTACAATCTTTTAG', 'MILSKRFTIF', 123123, 124324)
            possible_gene_start = start + putative_gene[0]
            possible_gene_end = start + putative_gene[1]

            possible_cds = SeqFeature(
                FeatureLocation(
                    possible_gene_start, possible_gene_end,
                    strand=putative_gene[2],
                ),
                type='CDS'
            )

            # Now we adjust our boundaries for the RBS that's required
            # There are only two cases, the rbs is upstream of it, or downstream
            if putative_gene[5] < possible_gene_start:
                possible_gene_start = putative_gene[5]
            else:
                possible_gene_end = putative_gene[6]

            possible_rbs = SeqFeature(
                FeatureLocation(
                    putative_gene[5], putative_gene[6],
                    strand=putative_gene[2],
                ),
                type='Shine_Dalgarno_sequence'
            )

            possible_gene = SeqFeature(
                FeatureLocation(
                    possible_gene_start, possible_gene_end,
                    strand=putative_gene[2],
                ),
                type='gene',
                qualifiers={
                    'note': ['Possible gene']
                }
            )
            possible_gene.sub_features = [possible_rbs, possible_cds]
            qc_features.append(possible_gene)

        better_results.append(result_obj + [len(putative_genes)])

    # Bad gaps are those with more than zero possible genes found
    bad = len([x for x in better_results if x[2] > 0])
    # Generally taking "good" here as every possible gap in the genome
    # Thus, good is TOTAL - gaps
    good = len(sorted_genes) + 1 - bad
    # and bad is just gaps
    return good, bad, better_results, qc_features


def phi(x):
    """Standard phi function used in calculation of normal distribution"""
    return math.exp(-1 * math.pi * x * x)


def norm(x, mean=0, sd=1):
    """
    Normal distribution. Given an x position, a mean, and a standard
    deviation, calculate the "y" value. Useful for score scaling

    Modified to multiply by SD. This means even at sd=5, norm(x, mean) where x = mean => 1, rather than 1/5.
    """
    return (1 / float(sd)) * phi(float(x - mean) / float(sd)) * sd


def coding_density(record, mean=92.5, sd=20):
    """
    Find coding density in the genome
    """
    feature_lengths = 0

    for gene_a in coding_genes(record.features):
        feature_lengths += sum([
            len(x) for x in
            genes(gene_a.sub_features, feature_type='CDS')
        ])

    avgFeatLen = float(feature_lengths) / float(len(record.seq))
    return int(norm(100 * avgFeatLen, mean=mean, sd=sd) * 100), int(100 * avgFeatLen)


def excessive_overlap(record, excess=15, excess_divergent=30):
    """
    Find excessive overlaps in the genome, where excessive is defined as 15
    bases for same strand, and 30 for divergent translation.

    Does a product of all the top-level features in the genome, and calculates
    gaps.
    """
    results = []
    bad = 0
    qc_features = []

    for (gene_a, gene_b) in itertools.combinations(coding_genes(record.features), 2):
        # Get the CDS from the subfeature list.
        # TODO: not recursive.
        cds_a = [x for x in genes(gene_a.sub_features, feature_type='CDS')]
        cds_b = [x for x in genes(gene_b.sub_features, feature_type='CDS')]

        if len(cds_a) == 0:
            log.warn("Gene missing subfeatures; %s", get_gff3_id(gene_a))
            continue

        if len(cds_b) == 0:
            log.warn("Gene missing subfeatures; %s", get_gff3_id(gene_b))
            continue

        cds_a = cds_a[0]
        cds_b = cds_b[0]

        # Set of locations that are included in the CDS of A and the
        # CDS of B
        cas = set(range(cds_a.location.start, cds_a.location.end))
        cbs = set(range(cds_b.location.start, cds_b.location.end))

        # Here we calculate the intersection between the two sets, and
        # if it's larger than our excessive size, we know that they're
        # overlapped
        ix = cas.intersection(cbs)

        if (cds_a.location.strand == cds_b.location.strand and len(ix) >= excess) or \
                (cds_a.location.strand != cds_b.location.strand and len(ix) >= excess_divergent):
            bad += float(len(ix)) / float(min(excess, excess_divergent))
            qc_features.append(gen_qc_feature(
                min(ix),
                max(ix),
                "Excessive Overlap",
                id_src=gene_a
            ))
            results.append((gene_a, gene_b, min(ix), max(ix)))

    # Good isn't accurate here. It's a triangle number and just ugly, but we
    # don't care enough to fix it.
    good = len(list(coding_genes(record.features)))
    good = int(good - bad)
    if good < 0:
        good = 0
    return good, int(bad), results, qc_features


def get_encouragement(score):
    """Some text telling the user how they did
    """
    for encouragement in ENCOURAGEMENT:
        if score > encouragement[0]:
            return encouragement[1]
    return ENCOURAGEMENT[-1][1]


def genome_overview(record):
    """Genome overview
    """
    data = {
        'genes': {
            'count': 0,
            'bases': 0,
            'density': 0,  # genes / kb
            'avg_len': [],
            'comp': {
                'A': 0,
                'C': 0,
                'G': 0,
                'T': 0,
            }
        },
        'overall': {
            'comp': {
                'A': record.seq.count('A'),
                'C': record.seq.count('C'),
                'G': record.seq.count('G'),
                'T': record.seq.count('T'),
            },
            'gc': 0,
        }
    }
    gene_features = list(coding_genes(record.features))
    data['genes']['count'] = len(gene_features)

    for feat in gene_features:
        data['genes']['comp']['A'] += feat.extract(record).seq.count('A')
        data['genes']['comp']['C'] += feat.extract(record).seq.count('C')
        data['genes']['comp']['T'] += feat.extract(record).seq.count('T')
        data['genes']['comp']['G'] += feat.extract(record).seq.count('G')
        data['genes']['bases'] += len(feat)
        data['genes']['avg_len'].append(len(feat))

    data['genes']['avg_len'] = float(sum(data['genes']['avg_len'])) / len(gene_features)
    data['overall']['gc'] = float(data['overall']['comp']['G'] + data['overall']['comp']['C']) / len(record.seq)
    return data


def find_morons(record):
    """Locate morons in the genome

    Don't even know why...

    TODO: remove? Idk.
    """
    results = []
    good = 0
    bad = 0

    gene_features = list(coding_genes(record.features))
    for i, gene in enumerate(gene_features):
        two_left = gene_features[i - 2:i]
        two_right = gene_features[i + 1:i + 1 + 2]
        strands = [x.strand for x in two_left] + [x.strand for x in two_right]
        anticon = [x for x in strands if x != gene.strand]

        if len(anticon) == 4:
            has_rbs = [x.type == "Shine_Dalgarno_sequence" for x in
                       gene.sub_features]
            if any(has_rbs):
                rbs = [x for x in gene.sub_features if x.type ==
                       "Shine_Dalgarno_sequence"][0]
                rbs_msg = str(rbs.extract(record.seq))
            else:
                rbs_msg = "No RBS Available"
            results.append((gene, two_left, two_right, rbs_msg))
            bad += 1
        else:
            good += 1
    return good, bad, results, []


def bad_gene_model(record):
    """Find features without product
    """
    results = []
    good = 0
    bad = 0
    qc_features = []

    for gene in coding_genes(record.features):
        exons = [x for x in genes(gene.sub_features, feature_type='exon') if len(x) > 10]
        CDSs = [x for x in genes(gene.sub_features, feature_type='CDS')]

        if len(exons) >= 1 and len(CDSs) >= 1:
            if len(exons) != len(CDSs):
                results.append((
                    get_gff3_id(gene),
                    None,
                    None,
                    'Mismatched number of exons and CDSs in gff3 representation',
                ))
                qc_features.append(gen_qc_feature(
                    gene.location.start, gene.location.end,
                    'Mismatched number of exons and CDSs in gff3 representation',
                    strand=gene.strand,
                    id_src=gene
                ))
                bad += 1
            else:
                for (exon, cds) in zip(exons, CDSs):
                    if len(exon) != len(cds):
                        results.append((
                            get_gff3_id(gene),
                            exon,
                            cds,
                            'CDS does not extend to full length of gene',
                        ))
                        qc_features.append(gen_qc_feature(
                            exon.location.start, exon.location.end,
                            'CDS does not extend to full length of gene',
                            strand=exon.strand,
                            id_src=gene
                        ))
                        bad += 1
                    else:
                        good += 1
        else:
            log.warn("Could not handle %s, %s", exons, CDSs)
            results.append((
                get_gff3_id(gene),
                None,
                None,
                '{0} exons, {1} CDSs'.format(len(exons), len(CDSs))
            ))

    return good, len(results) + bad, results, qc_features


def weird_starts(record):
    """Find features without product
    """
    good = 0
    bad = 0
    qc_features = []
    results = []

    overall = {}
    for gene in coding_genes(record.features):
        seq = [x for x in genes(gene.sub_features, feature_type='CDS')]
        if len(seq) == 0:
            log.warn("No CDS for gene %s", get_gff3_id(gene))
            continue
        else:
            seq = seq[0]

        seq_str = str(seq.extract(record.seq))
        start_codon = seq_str[0:3]
        stop_codon = seq_str[-3]
        seq.__start = start_codon
        seq.__stop = stop_codon
        if start_codon not in overall:
            overall[start_codon] = 1
        else:
            overall[start_codon] += 1

        if start_codon not in ('ATG', 'TTG', 'GTG'):
            log.warn("Weird start codon (%s) on %s", start_codon, get_gff3_id(gene))
            seq.__error = 'Unusual start codon %s' % start_codon

            s = 0
            e = 0
            if seq.strand > 0:
                s = seq.location.start
                e = seq.location.start + 3
            else:
                s = seq.location.end
                e = seq.location.end - 3

            results.append(seq)

            qc_features.append(gen_qc_feature(
                s, e,
                'Weird start codon',
                strand=seq.strand,
                id_src=gene
            ))
            bad += 1
        else:
            good += 1

    return good, bad, results, qc_features, overall


def missing_genes(record):
    """Find features without product
    """
    results = []
    good = 0
    bad = 0
    qc_features = []

    for gene in coding_genes(record.features):
        if gene.qualifiers.get('cpt_source', [None])[0] == 'CPT_GENE_MODEL_CORRECTION':
            results.append(gene)
            bad += 1
        else:
            good += 1

    return good, bad, results, qc_features


def gene_model_correction_issues(record):
    """Find features that have issues from the gene model correction step.
    These have qualifiers beginning with CPT_GMS
    """
    results = []
    good = 0
    bad = 0
    qc_features = []

    # For each gene
    for gene in coding_genes(record.features):
        # Get the list of child CDSs
        cdss = [x for x in genes(gene.sub_features, feature_type='CDS')]
        # And our matching qualifiers
        gene_data = [(k, v) for (k, v) in gene.qualifiers.items() if k == 'cpt_gmc']
        # If there are problems with ONLY the parent, let's complain
        local_results = []
        local_qc_features = []
        for x in gene_data:
            if 'Missing Locus Tag' in x[1]:
                # Missing locus tag is an either or thing, if it hits here
                # there shouldn't be anything else wrong with it.

                # Obviously missing so we remove it
                gene.qualifiers['locus_tag'] = [""]
                # Translation from bp_genbank2gff3.py
                cdss[0].qualifiers['locus_tag'] = cdss[0].qualifiers['Name']
                # Append our results
                local_results.append((
                    gene,
                    cdss[0],
                    'Gene is missing a locus_tag'
                ))
                local_qc_features.append(gen_qc_feature(
                    gene.location.start,
                    gene.location.end,
                    'Gene is missing a locus_tag',
                    strand=gene.strand
                ))

        # We need to alert on any child issues as well.
        for cds in cdss:
            cds_data = [(k, v[0]) for (k, v) in cds.qualifiers.items() if k == 'cpt_gmc']
            if len(gene_data) == 0 and len(cds_data) == 0:
                # Alles gut
                pass
            else:
                for _, problem in cds_data:
                    if problem == 'BOTH Missing Locus Tag':
                        gene.qualifiers['locus_tag'] = ['']
                        cds.qualifiers['locus_tag'] = ['']
                        local_results.append((
                            gene, cds,
                            'Both gene and CDS are missing locus tags'
                        ))
                        local_qc_features.append(gen_qc_feature(
                            cds.location.start,
                            cds.location.end,
                            'CDS is missing a locus_tag',
                            strand=cds.strand
                        ))
                        local_qc_features.append(gen_qc_feature(
                            gene.location.start,
                            gene.location.end,
                            'Gene is missing a locus_tag',
                            strand=gene.strand
                        ))
                    elif problem == 'Different locus tag from associated gene.':
                        gene.qualifiers['locus_tag'] = gene.qualifiers['Name']
                        cds.qualifiers['locus_tag'] = cds.qualifiers['cpt_gmc_locus']
                        local_results.append((
                            gene, cds,
                            'Gene and CDS have differing locus tags',
                        ))
                        local_qc_features.append(gen_qc_feature(
                            gene.location.start,
                            gene.location.end,
                            'Gene and CDS have differing locus tags',
                            strand=gene.strand
                        ))
                    elif problem == 'Missing Locus Tag':
                        # Copy this over
                        gene.qualifiers['locus_tag'] = gene.qualifiers['Name']
                        # This one is missing
                        cds.qualifiers['locus_tag'] = ['']
                        local_results.append((
                            gene, cds,
                            'CDS is missing a locus_tag',
                        ))
                        local_qc_features.append(gen_qc_feature(
                            cds.location.start,
                            cds.location.end,
                            'CDS is missing a locus_tag',
                            strand=cds.strand
                        ))
                    else:
                        log.warn("Cannot handle %s", problem)

        if len(local_results) > 0:
            bad += 1
        else:
            good += 1

        qc_features.extend(local_qc_features)
        results.extend(local_results)
    return good, bad, results, qc_features


def missing_tags(record):
    """Find features without product
    """
    results = []
    good = 0
    bad = 0
    qc_features = []

    for gene in coding_genes(record.features):
        cds = [x for x in genes(gene.sub_features, feature_type='CDS')]
        if len(cds) == 0:
            log.warn("Gene missing CDS subfeature %s", get_gff3_id(gene))
            continue

        cds = cds[0]

        if 'product' not in cds.qualifiers:
            log.info("Missing product tag on %s", get_gff3_id(gene))
            qc_features.append(gen_qc_feature(
                cds.location.start,
                cds.location.end,
                'Missing product tag',
                strand=cds.strand
            ))
            results.append(cds)
            bad += 1
        else:
            good += 1

    return good, bad, results, qc_features


def evaluate_and_report(annotations, genome, gff3=None,
                        tbl=None, sd_min=5, sd_max=15, min_gene_length=30,
                        excessive_gap_dist=50, excessive_gap_divergent_dist=200,
                        excessive_overlap_dist=25, excessive_overlap_divergent_dist=50,
                        reportTemplateName='phage_annotation_validator.html'):
    """
    Generate our HTML evaluation of the genome
    """
    # Get features from GFF file
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    # Get the first GFF3 record
    # TODO: support multiple GFF3 files.
    record = list(GFF.parse(annotations, base_dict=seq_dict))[0]

    gff3_qc_record = SeqRecord(record.id, id=record.id)
    gff3_qc_record.features = []
    gff3_qc_features = []



    log.info("Locating missing RBSs")
    # mb_any = "did they annotate ANY rbss? if so, take off from score."
    mb_good, mb_bad, mb_results, mb_annotations, mb_any = missing_rbs(
        record,
        lookahead_min=sd_min,
        lookahead_max=sd_max
    )
    gff3_qc_features += mb_annotations

    log.info("Locating excessive gaps")
    eg_good, eg_bad, eg_results, eg_annotations = excessive_gap(
        record,
        excess=excessive_gap_dist,
        excess_divergent=excessive_gap_divergent_dist,
        min_gene=min_gene_length,
        slop=excessive_overlap_dist,
        lookahead_min=sd_min,
        lookahead_max=sd_max
    )
    gff3_qc_features += eg_annotations

    log.info("Locating excessive overlaps")
    eo_good, eo_bad, eo_results, eo_annotations = excessive_overlap(
        record,
        excess=excessive_overlap_dist,
        excess_divergent=excessive_overlap_divergent_dist
    )
    gff3_qc_features += eo_annotations

    log.info("Locating morons")
    mo_good, mo_bad, mo_results, mo_annotations = find_morons(record)
    gff3_qc_features += mo_annotations

    log.info("Locating missing tags")
    mt_good, mt_bad, mt_results, mt_annotations = missing_tags(record)
    gff3_qc_features += mt_annotations

    log.info("Locating missing gene features")
    mg_good, mg_bad, mg_results, mg_annotations = missing_genes(record)
    gff3_qc_features += mg_annotations

    log.info("Determining coding density")
    cd, cd_real = coding_density(record)

    log.info("Locating weird starts")
    ws_good, ws_bad, ws_results, ws_annotations, ws_overall = weird_starts(record)
    gff3_qc_features += ws_annotations

    log.info("Locating bad gene models")
    gm_good, gm_bad, gm_results, gm_annotations = bad_gene_model(record)
    if gm_good + gm_bad == 0:
        gm_bad = 1

    log.info("Locating more bad gene models")
    gmc_good, gmc_bad, gmc_results, gmc_annotations = gene_model_correction_issues(record)
    if gmc_good + gmc_bad == 0:
        gmc_bad = 1

    good_scores = [eg_good, eo_good, mt_good, ws_good, gm_good, gmc_good]
    bad_scores = [eg_bad, eo_bad, mt_bad, ws_bad, gm_bad, gmc_bad]

    # Only if they tried to annotate RBSs do we consider them.
    if mb_any:
        good_scores.append(mb_good)
        bad_scores.append(mb_bad)
    subscores = []

    for (g, b) in zip(good_scores, bad_scores):
        if g + b == 0:
            s = 0
        else:
            s = int(100 * float(g) / (float(b) + float(g)))
        subscores.append(s)
    subscores.append(cd)

    score = int(float(sum(subscores)) / float(len(subscores)))

    # This is data that will go into our HTML template
    kwargs = {
        'upstream_min': sd_min,
        'upstream_max': sd_max,
        'record_name': record.id,
        'record_nice_name': nice_name(record),
        'params': {
            'sd_min': sd_min,
            'sd_max': sd_max,
            'min_gene_length': min_gene_length,
            'excessive_gap_dist': excessive_gap_dist,
            'excessive_gap_divergent_dist': excessive_gap_divergent_dist,
            'excessive_overlap_dist': excessive_overlap_dist,
            'excessive_overlap_divergent_dist': excessive_overlap_divergent_dist,
        },
        'score': score,
        'encouragement': get_encouragement(score),
        'genome_overview' : genome_overview(record),

        'rbss_annotated': mb_any,
        'missing_rbs': mb_results,
        'missing_rbs_good': mb_good,
        'missing_rbs_bad': mb_bad,
        'missing_rbs_score': 0 if mb_good + mb_bad == 0 else (100 * mb_good / (mb_good + mb_bad)),

        'excessive_gap': eg_results,
        'excessive_gap_good': eg_good,
        'excessive_gap_bad': eg_bad,
        'excessive_gap_score': 0 if eo_good + eo_bad == 0 else (100 * eo_good / (eo_good + eo_bad)),

        'excessive_overlap': eo_results,
        'excessive_overlap_good': eo_good,
        'excessive_overlap_bad': eo_bad,
        'excessive_overlap_score': 0 if eo_good + eo_bad == 0 else (100 * eo_good / (eo_good + eo_bad)),

        'morons': mo_results,
        'morons_good': mo_good,
        'morons_bad': mo_bad,
        'morons_score': 0 if mo_good + mo_bad == 0 else (100 * mo_good / (mo_good + mo_bad)),

        'missing_tags': mt_results,
        'missing_tags_good': mt_good,
        'missing_tags_bad': mt_bad,
        'missing_tags_score': 0 if mt_good + mt_bad == 0 else (100 * mt_good / (mt_good + mt_bad)),

        'missing_genes': mg_results,
        'missing_genes_good': mg_good,
        'missing_genes_bad': mg_bad,
        'missing_genes_score': 0 if mg_good + mg_bad == 0 else (100 * mg_good / (mg_good + mg_bad)),

        'weird_starts': ws_results,
        'weird_starts_good': ws_good,
        'weird_starts_bad': ws_bad,
        'weird_starts_overall': ws_overall,
        'weird_starts_overall_sorted_keys': sorted(ws_overall, reverse=True, key=lambda x: ws_overall[x]),
        'weird_starts_score': 0 if ws_good + ws_bad == 0 else (100 * ws_good / (ws_good + ws_bad)),

        'gene_model': gm_results,
        'gene_model_good': gm_good,
        'gene_model_bad': gm_bad,
        'gene_model_score': 0 if gm_good + gm_bad == 0 else (100 * gm_good / (gm_good + gm_bad)),

        'gene_model_correction': gmc_results,
        'gene_model_correction_good': gmc_good,
        'gene_model_correction_bad': gmc_bad,
        'gene_model_correction_score': 0 if gmc_good + gmc_bad == 0 else (100 * gmc_good / (gmc_good + gmc_bad)),

        'coding_density': cd,
        'coding_density_real': cd_real,
        'coding_density_score': cd,
    }

    with open(tbl, 'w') as handle:
        kw_subset = {}
        for key in kwargs:
            if key in ('score', 'record_name') or '_good' in key or '_bad' in key or '_overall' in key:
                kw_subset[key] = kwargs[key]
        json.dump(kw_subset, handle)

    with open(gff3, 'w') as handle:
        gff3_qc_record.features = gff3_qc_features
        gff3_qc_record.annotations = {}
        GFF.write([gff3_qc_record], handle)

    def nice_strand(direction):
        if direction > 0:
            return '→'.decode('utf-8')
        else:
            return '←'.decode('utf-8')

    def nice_strand_tex(direction):
        if direction > 0:
            return '$\\rightarrow$'
        else:
            return '$\\leftarrow$'

    def texify(data):
        return data.replace('_', '\\_').replace('$', '\\$')

    def length(data):
        return len(data)

    env = Environment(loader=FileSystemLoader(SCRIPT_PATH), trim_blocks=True, lstrip_blocks=True)
    env.filters.update({
        'nice_id': get_gff3_id,
        'nice_strand': nice_strand,
        'nice_strand_tex': nice_strand_tex,
        'texify': texify,
        'length': length,
    })
    tpl = env.get_template(reportTemplateName)
    return tpl.render(**kwargs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('annotations', type=argparse.FileType("r"), help='Parent GFF3 annotations')
    parser.add_argument('genome', type=argparse.FileType("r"), help='Genome Sequence')
    parser.add_argument('--gff3', type=str, help='GFF3 Annotations', default='qc_annotations.gff3')
    parser.add_argument('--tbl', type=str, help='Table for noninteractive parsing', default='qc_results.json')

    parser.add_argument('--sd_min', type=int, help='Minimum distance from gene start for an SD to be', default=5)
    parser.add_argument('--sd_max', type=int, help='Maximum distance from gene start for an SD to be', default=15)

    parser.add_argument('--min_gene_length', type=int, help='Minimum length for a putative gene call (AAs)', default=30)

    parser.add_argument('--excessive_overlap_dist', type=int, help='Excessive overlap for genes in same direction', default=25)
    parser.add_argument('--excessive_overlap_divergent_dist', type=int, help='Excessive overlap for genes in diff directions', default=50)

    parser.add_argument('--excessive_gap_dist', type=int, help='Maximum distance between two genes', default=40)
    parser.add_argument('--excessive_gap_divergent_dist', type=int, help='Maximum distance between two divergent genes', default=200)

    parser.add_argument('--reportTemplateName', help='Report template file name', default='phageqc_report_full.html')

    args = parser.parse_args()

    print evaluate_and_report(**vars(args))
    # evaluate_and_report(**vars(args))
