#!/usr/bin/env python
import sys
import math
import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from jinja2 import Template
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

REPORT_TEMPLATE = Template(open('phage_annotation_validator.html').read())
dashcss = open('dashboard.css').read()

__author__ = "Eric Rasche"
__version__ = "0.4.0"
__maintainer__ = "Eric Rasche"
__email__ = "esr@tamu.edu"

ENCOURAGEMENT = (
    (100, 'Perfection itself!'),
    (90, 'Not too bad, a few minor things to fix...'),
    (50, 'Issues detected! </p><p class="text-muted">Have you heard of the <a href="https://cpt.tamu.edu">CPT</a>\'s Automated Phage Annotation Pipeline?'),
    (0, '<b>MAJOR</b> issues detected! Please strongly consider using the <a href="https://cpt.tamu.edu">CPT</a>\'s Automated Phage Annotation Pipeline'),
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
    results = []
    good = 0
    bad = 0

    for gene in genes(record.features):
        has_rbs = [x.type == "Shine_Dalgarno_sequence" for x in gene.sub_features]
        if not any(has_rbs):
            # Shinefind copy pasta to find upstream seq
            if gene.strand > 0:
                start = gene.location.start - lookahead_max
                end = gene.location.start - lookahead_min
            else:
                start = gene.location.end + lookahead_min
                end = gene.location.end + lookahead_max
            (start, end) = __ensure_location_in_bounds(start=start, end=end,
                                                       parent_length=record.__len__)
            tmp = SeqFeature(FeatureLocation(start, end, strand=gene.strand),
                             type='domain')
            seq = str(tmp.extract(record.seq))
            gene.__upstream = seq
            gene.__message = "No RBS"

            bad += 1
            results.append(gene)
        else:
            rbs = [x for x in gene.sub_features if x.type == "Shine_Dalgarno_sequence"][0]
            cds = [x for x in gene.sub_features if x.type == "CDS"][0]

            if gene.strand > 0:
                distance = cds.location.start - rbs.location.end
            else:
                distance = rbs.location.start - cds.location.end

            if distance > lookahead_max:
                gene.__message = "RBS too far away (%s nt)" % distance
                bad += 1
                results.append(gene)
            else:
                good += 1

    return good, bad, results


def excessive_gap(record, excess=10):
    results = []
    good = 0
    bad = 0

    annotated_regions = {}
    for gene in genes(record.features):
        for i in range(gene.location.start, gene.location.end):
            annotated_regions[i] = True

    annotated = True
    unannotated_count = 0
    region_start = 0
    for i in range(1, len(record.seq)):
        # If it's not annotated
        if i not in annotated_regions:
            # Bump the count
            unannotated_count += 1
            # Moving from annotated to unannotated
            if annotated:
                region_start = i
            # Ensure we mark that we're in an unannotated region
            annotated = False
        else:
            # IF we're switching to an annotated state from an unannotated
            # state...
            if not annotated:
                if unannotated_count > excess:
                    results.append((region_start, i))
                unannotated_count = 0

            annotated = True

    if not annotated:
        results.append((region_start, i))

    # Generally taking "good" here as every possible gap in the genome
    good = len(list(genes(record.features))) + 1 - len(results)
    # Thus, good is TOTAL - gaps
    bad = len(results)
    # and bad is just gaps
    return good, bad, results


def genes(feature_list):
    for x in feature_list:
        if x.type == 'gene':
            yield x


def excessive_overlap(record, excessive=15):
    results = []
    bad = 0

    for gene_a in genes(record.features):
        for gene_b in genes(record.features):
            if gene_a != gene_b and gene_a.location.start <= gene_b.location.start:
                cds_a = [x for x in gene_a.sub_features if x.type == "CDS"][0]
                cds_b = [x for x in gene_b.sub_features if x.type == "CDS"][0]

                cas = set(range(cds_a.location.start, cds_a.location.end))
                cbs = set(range(cds_b.location.start, cds_b.location.end))

                ix = cas.intersection(cbs)
                if len(ix) >= excessive:
                    bad += 1
                    results.append((gene_a, gene_b, min(ix), max(ix)))

    # Good isn't accurate here. It's a triangle number and just ugly, but we
    # don't care enough to fix it.
    return len(list(genes(record.features))), bad, results


def get_encouragement(score):
    for encouragement in ENCOURAGEMENT:
        if score > encouragement[0]:
            return encouragement[1]
    return ENCOURAGEMENT[-1][1]


def find_morons(record):
    results = []
    good = 0
    bad = 0

    gene_features = list(genes(record.features))
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
    return good, bad, results


def evaluate_and_report(annotations, genome):
    # Get features from GFF file
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    record = list(GFF.parse(annotations, base_dict=seq_dict))[0]
    upstream_min = 5
    upstream_max = 15

    mb_good, mb_bad, mb_results = missing_rbs(record,
                                              lookahead_min=upstream_min,
                                              lookahead_max=upstream_max)
    eg_good, eg_bad, eg_results = excessive_gap(record, excess=3 * upstream_max)
    eo_good, eo_bad, eo_results = excessive_overlap(record, excessive=15)
    mo_good, mo_bad, mo_results = find_morons(record)

    score_good = float(sum((mb_good, eg_good, eo_good)))
    score_bad = float(sum((mb_bad, eg_bad, eo_bad)))

    score = int(100 * score_good / (score_bad + score_good))

    kwargs = {
        'dash': dashcss,
        'upstream_min': upstream_min,
        'upstream_max': upstream_max,
        'record_name': record.id,

        'score': score,
        'encouragement': get_encouragement(score),

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
    }

    return REPORT_TEMPLATE.render(**kwargs)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('annotations', type=file, help='Parent GFF3 annotations')
    parser.add_argument('genome', type=file, help='Genome Sequence')
    args = parser.parse_args()

    print evaluate_and_report(**vars(args))
