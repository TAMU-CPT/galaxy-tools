#!/usr/bin/env python
import os
import json
import math
import argparse
import itertools
from gff3 import feature_lambda, feature_test_type, feature_test_quals
from shinefind import NaiveSDCaller
from BCBio import GFF
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import reverse_complement, translate
from Bio.SeqFeature import SeqFeature, FeatureLocation
from jinja2 import Template
import re
import logging
logging.basicConfig(level=logging.WARN)
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


def gen_qc_feature(start, end, message, strand=0):
    return SeqFeature(
        FeatureLocation(start, end, strand=strand),
        qualifiers={
            'note': [message]
        }
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
        rbs_rbs = list(feature_lambda(gene.sub_features, feature_test_type, {'type': 'RBS'}, subfeatures=False))
        rbs_sds = list(feature_lambda(gene.sub_features,
                                      feature_test_type, {'type': 'Shine_Dalgarno_sequence'}, subfeatures=False))
        regulatory_elements = list(feature_lambda(gene.sub_features,
                                                  feature_test_type, {'type': 'regulatory'}, subfeatures=False))
        rbs_regulatory = list(feature_lambda(regulatory_elements,
                                             feature_test_quals, {'regulatory_class': ['ribosome_binding_site']},
                                             subfeatures=False))

        rbss = rbs_rbs + rbs_sds + rbs_regulatory
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


            qc_features.append(gen_qc_feature(start, end, 'Missing RBS', strand=gene.strand))

            bad += 1
            results.append(gene)
        else:
            if len(rbss) > 1:
                log.warn("%s RBSs found for gene %s", rbss[0].id, gene.id)
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
                    strand=gene.strand))

                bad += 1
                results.append(gene)
            else:
                good += 1

    return good, bad, results, qc_features, any_rbss

# modified from get_orfs_or_cdss.py
# -----------------------------------------------------------

# get stop codons
table_obj = CodonTable.ambiguous_generic_by_id[11]

starts = sorted(table_obj.start_codons)
re_starts = re.compile("|".join(starts))

stops = sorted(table_obj.stop_codons)
re_stops = re.compile("|".join(stops))


def start_chop_and_trans(s, strict=True):
    """Returns offset, trimmed nuc, protein."""
    if strict:
        assert s[-3:] in stops, s
    assert len(s) % 3 == 0
    for match in re_starts.finditer(s):
        # Must check the start is in frame
        start = match.start()
        if start % 3 == 0:
            n = s[start:]
            assert len(n) % 3 == 0, "%s is len %i" % (n, len(n))
            if strict:
                t = translate(n, 11, cds=True)
            else:
                # Use when missing stop codon,
                t = "M" + translate(n[3:], 11, to_stop=True)
            return start, n, t
    return None, None, None


def break_up_frame(s, min_len=10):
    """Returns offset, nuc, protein."""
    start = 0
    for match in re_stops.finditer(s):
        index = match.start() + 3
        if index % 3 != 0:
            continue
        n = s[start:index]

        offset, n, t = start_chop_and_trans(n)
        if n and len(t) >= min_len:
            yield start + offset, n, t
        start = index


def putative_genes_in_sequence(nuc_seq, min_len=10):
    """Returns start, end, strand, nucleotides, protein.
    Co-ordinates are Python style zero-based.
    """
    nuc_seq = nuc_seq.upper()
    # TODO - Refactor to use a generator function (in start order)
    # rather than making a list and sorting?
    answer = []
    full_len = len(nuc_seq)

    for frame in range(0, 3):
        for offset, n, t in break_up_frame(nuc_seq[frame:], min_len=min_len):
            start = frame + offset  # zero based
            answer.append((start, start + len(n), +1, n, t))

    rc = reverse_complement(nuc_seq)
    for frame in range(0, 3):
        for offset, n, t in break_up_frame(rc[frame:], min_len=min_len):
            start = full_len - frame - offset  # zero based
            answer.append((start, start - len(n), -1, n, t))
    answer.sort()
    return answer


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


def excessive_gap(record, excess=10, min_gene=30, slop=30, lookahead_min=5, lookahead_max=15):
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

    current_gene = [
        int(sorted_genes[0].location.start),
        int(sorted_genes[0].location.end)
    ]
    for gene in sorted_genes[1:]:
        # If the gene's start is contiguous to the "current_gene", then we
        # extend current_gene
        if gene.location.start <= current_gene[1] + excess:
            current_gene[1] = int(gene.location.end)
        else:
            # If it's discontiguous, we append the region and clear.
            contiguous_regions.append(current_gene)
            current_gene = [int(gene.location.start), int(gene.location.end)]

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

        if b[0] > a[1] + excess:
            results.append((a[1], b[0]))

    better_results = []
    qc_features = []
    for (start, end) in results:
        f = gen_qc_feature(start, end, 'Excessive gap, %s bases' % abs(end-start))
        qc_features.append(f)
        putative_genes = putative_genes_in_sequence(str(record[start - slop:end + slop].seq), min_len=min_gene)
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

        better_results.append((start, end, len(putative_genes)))

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


def coding_genes(feature_list):
    for x in feature_lambda(feature_list, feature_test_type, {'type': 'gene'}, subfeatures=True):
        if len(list(feature_lambda(x.sub_features, feature_test_type, {'type': 'CDS'}, subfeatures=False))) > 0:
            yield x


def genes(feature_list, feature_type='gene'):
    """
    Simple filter to extract gene features from the feature set.
    """

    for x in feature_lambda(feature_list, feature_test_type,
                            {'type': feature_type},
                            subfeatures=True):
        yield x


def excessive_overlap(record, excessive=15):
    """
    Find excessive overlaps in the genome, where excessive is defined as 15
    bases.

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
            log.warn("Gene missing subfeatures; %s", gene_a.id)
            continue

        if len(cds_b) == 0:
            log.warn("Gene missing subfeatures; %s", gene_b.id)
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
        if len(ix) >= excessive:
            bad += float(len(ix)) / float(excessive)
            qc_features.append(gen_qc_feature(
                min(ix),
                max(ix),
                "Excessive Overlap")
            )
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


def weird_starts(record):
    """Find features without product
    """
    results = []
    good = 0
    bad = 0
    qc_features = []

    overall = {}
    for gene in coding_genes(record.features):
        seq = [x for x in genes(gene.sub_features, feature_type='CDS')]
        if len(seq) == 0:
            log.warn("No CDS for gene %s", gene.id)
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
            log.warn("Weird start codon (%s) on %s", start_codon, seq.id)
            seq.__error = 'Unusual start codon %s' % start_codon

            s = 0
            e = 0
            if seq.strand > 0:
                s = seq.location.start
                e = seq.location.start + 3
            else:
                s = seq.location.end
                e = seq.location.end - 3

            qc_features.append(gen_qc_feature(
                s, e,
                'Weird start codon',
                strand=seq.strand
            ))
            results.append(seq)
            bad += 1
        else:
            good += 1

    return good, bad, results, qc_features, overall


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
            log.warn("Gene missing CDS subfeature %s", gene.id)
            continue

        cds = cds[0]

        if 'product' not in cds.qualifiers:
            log.warn("Missing product tag on %s", cds.id)
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


def evaluate_and_report(annotations, genome, user_email, gff3=None, tbl=None, sd_min=5,
                        sd_max=15, gap_dist=45, overlap_dist=15,
                        min_gene_length=30, reportTemplateName='phage_annotation_validator.html'):
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
    #mb_any = "did they annotate ANY rbss? if so, take off from score."
    mb_good, mb_bad, mb_results, mb_annotations, mb_any = missing_rbs(
        record,
        lookahead_min=sd_min,
        lookahead_max=sd_max
    )
    gff3_qc_features += mb_annotations

    log.info("Locating excessive gaps")
    eg_good, eg_bad, eg_results, eg_annotations = excessive_gap(
        record, excess=gap_dist, min_gene=min_gene_length,
        slop=overlap_dist, lookahead_min=sd_min, lookahead_max=sd_max)
    gff3_qc_features += eg_annotations

    log.info("Locating excessive overlaps")
    eo_good, eo_bad, eo_results, eo_annotations = excessive_overlap(record, excessive=overlap_dist)
    gff3_qc_features += eo_annotations

    log.info("Locating morons")
    mo_good, mo_bad, mo_results, mo_annotations = find_morons(record)
    gff3_qc_features += mo_annotations

    log.info("Locating missing tags")
    mt_good, mt_bad, mt_results, mt_annotations = missing_tags(record)
    gff3_qc_features += mt_annotations

    log.info("Determining coding density")
    cd, cd_real = coding_density(record)

    log.info("Locating weird starts")
    ws_good, ws_bad, ws_results, ws_annotations, ws_overall = weird_starts(record)
    gff3_qc_features += ws_annotations


    good_scores = [eg_good, eo_good, mt_good, ws_good]
    bad_scores = [eg_bad, eo_bad, mt_bad, ws_bad]

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

    from guanine import GuanineClient
    g = GuanineClient()
    g.submit(user_email, 'WfB-PhageQC', float(score) / 100)

    # This is data that will go into our HTML template
    kwargs = {
        'upstream_min': sd_min,
        'upstream_max': sd_max,
        'record_name': record.id,

        'score': score,
        'encouragement': get_encouragement(score),

        'rbss_annotated': mb_any,
        'missing_rbs': mb_results,
        'missing_rbs_good': mb_good,
        'missing_rbs_bad': mb_bad,

        'excessive_gap': eg_results,
        'excessive_gap_good': eg_good,
        'excessive_gap_bad': eg_bad,

        'excessive_overlap': eo_results,
        'excessive_overlap_good': eo_good,
        'excessive_overlap_bad': eo_bad,

        'morons': mo_results,
        'morons_good': mo_good,
        'morons_bad': mo_bad,

        'missing_tags': mt_results,
        'missing_tags_good': mt_good,
        'missing_tags_bad': mt_bad,

        'weird_starts': ws_results,
        'weird_starts_good': ws_good,
        'weird_starts_bad': ws_bad,
        'weird_starts_overall': ws_overall,
        'weird_starts_overall_sorted_keys': sorted(ws_overall, reverse=True, key=lambda x: ws_overall[x]),

        'coding_density': cd,
        'coding_density_real': cd_real,
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

    report_template = Template(open(os.path.join(SCRIPT_PATH, reportTemplateName), 'r').read())
    return report_template.render(**kwargs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('annotations', type=file, help='Parent GFF3 annotations')
    parser.add_argument('genome', type=file, help='Genome Sequence')
    parser.add_argument('--gff3', type=str, help='GFF3 Annotations', default='qc_annotations.gff3')
    parser.add_argument('--tbl', type=str, help='Table for noninteractive parsing', default='qc_results.json')

    parser.add_argument('--sd_min', type=int, help='Minimum distance from gene start for an SD to be', default=5)
    parser.add_argument('--sd_max', type=int, help='Maximum distance from gene start for an SD to be', default=15)

    parser.add_argument('--gap_dist', type=int, help='Maximum distance from gene start for an SD to be', default=30)
    parser.add_argument('--overlap_dist', type=int, help='Maximum distance from gene start for an SD to be', default=30)

    parser.add_argument('--min_gene_length', type=int, help='Minimum length for a putative gene call (AAs)', default=30)

    parser.add_argument('--reportTemplateName', help='Report template file name')
    parser.add_argument('--user_email')

    args = parser.parse_args()

    print evaluate_and_report(**vars(args))
