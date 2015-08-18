#!/usr/bin/env python
import os
import copy
import argparse
import itertools
from BCBio import GFF
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.Seq import reverse_complement, translate
from Bio.SeqFeature import SeqFeature, FeatureLocation
from jinja2 import Template
import logging
import re
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# Path to script, required because of Galaxy.
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
# Path to the HTML template for the report
REPORT_TEMPLATE = Template(open(os.path.join(SCRIPT_PATH, 'phage_annotation_validator.html'), 'r').read())

__author__ = "Eric Rasche"
__version__ = "0.4.0"
__maintainer__ = "Eric Rasche"
__email__ = "esr@tamu.edu"

ENCOURAGEMENT = (
    (100, 'Perfection itself!'),
    (90, 'Not too bad, a few minor things to fix...'),
    (70, 'Some issues to address'),
    (50, 'Issues detected! </p><p class="text-muted">Have you heard of the <a href="https://cpt.tamu.edu">CPT</a>\'s Automated Phage Annotation Pipeline?'),
    (0, '<b>MAJOR</b> issues detected! Please strongly consider using the <a href="https://cpt.tamu.edu">CPT</a>\'s Automated Phage Annotation Pipeline'),
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

    for gene in genes(record.features):
        # Check if there are RBSs, TODO: make this recursive. Each feature in
        # gene.sub_features can also have sub_features.
        has_rbs = [x.type == "Shine_Dalgarno_sequence" for x in gene.sub_features]
        # No RBS found
        if not any(has_rbs):
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
            gene.__upstream = seq
            gene.__message = "No RBS"

            qc_features.append(gen_qc_feature(start, end, 'Missing RBS', strand=gene.strand))

            bad += 1
            results.append(gene)
        else:
            # get first RBS/CDS
            rbs = [x for x in gene.sub_features if x.type == "Shine_Dalgarno_sequence"][0]
            cds = [x for x in gene.sub_features if x.type == "CDS"][0]

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

    return good, bad, results, qc_features

# modified from get_orfs_or_cdss.py
#-----------------------------------------------------------

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
        #Must check the start is in frame
        start = match.start()
        if start % 3 == 0:
            n = s[start:]
            assert len(n) % 3 == 0, "%s is len %i" % (n, len(n))
            if strict:
                t = translate(n, 11, cds=True)
            else:
                #Use when missing stop codon,
                t = "M" + translate(n[3:], 11, to_stop=True)
            return start, n, t
    return None, None, None


def break_up_frame(s):
    """Returns offset, nuc, protein."""
    start = 0
    for match in re_stops.finditer(s):
        index = match.start() + 3
        if index % 3 != 0:
            continue
        n = s[start:index]

        offset, n, t = start_chop_and_trans(n)
        if n and len(t) >= 10:
            yield start + offset, n, t
        start = index


def get_peptides(nuc_seq):
    """Returns start, end, strand, nucleotides, protein.
    Co-ordinates are Python style zero-based.
    """
    #TODO - Refactor to use a generator function (in start order)
    #rather than making a list and sorting?
    answer = []
    full_len = len(nuc_seq)

    for frame in range(0,3):
        for offset, n, t in break_up_frame(nuc_seq[frame:]):
            start = frame + offset #zero based
            answer.append((start, start + len(n), +1, n, t))

    rc = reverse_complement(nuc_seq)
    for frame in range(0,3) :
        for offset, n, t in break_up_frame(rc[frame:]):
            start = full_len - frame - offset #zero based
            answer.append((start - len(n), start, -1, n ,t))
    answer.sort()
    return answer


def putative_genes_in_sequence(sequence):
    return len(get_peptides(sequence.upper())) > 0


def excessive_gap(record, excess=10):
    """
    Identify excessive gaps between gene features.

    Default "excessive" gap size is 10, but that should likely be larger.
    """
    results = []
    good = 0
    bad = 0


    contiguous_regions = []
    sorted_genes = sorted(genes(record.features), key=lambda feature: feature.location.start)
    current_gene = [
        int(sorted_genes[0].location.start),
        int(sorted_genes[0].location.end)
    ]
    for gene in sorted_genes[1:]:
        if gene.location.start <= current_gene[1] + excess:
            current_gene[1] = int(gene.location.end)
        else:
            contiguous_regions.append(current_gene)
            current_gene = [int(gene.location.start), int(gene.location.end)]

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

    qc_features = [
        gen_qc_feature(
            start,
            end,
            'Excessive gap, %s bases' % abs(end-start))
        for (start, end) in results
    ]

    results = [(start, end, putative_genes_in_sequence(str(record[start:end].seq))) for (start, end) in results]
    # Bad gaps are those with more than zero possible genes found
    bad = len([x for x in results if x[2] > 0])
    # Generally taking "good" here as every possible gap in the genome
    # Thus, good is TOTAL - gaps
    good = len(list(genes(record.features))) + 1 - bad
    # and bad is just gaps
    return good, bad, results, qc_features


def genes(feature_list):
    """
    Simple filter to extract gene features from the feature set.

    TODO: not recursive.
    """
    for x in feature_list:
        if x.type == 'gene':
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

    for (gene_a, gene_b) in itertools.combinations(genes(record.features), 2):
        # Get the CDS from the subfeature list.
        # TODO: not recursive.
        cds_a = [x for x in gene_a.sub_features if x.type == "CDS"][0]
        cds_b = [x for x in gene_b.sub_features if x.type == "CDS"][0]

        # Set of locations that are included in the CDS of A and the
        # CDS of B
        cas = set(range(cds_a.location.start, cds_a.location.end))
        cbs = set(range(cds_b.location.start, cds_b.location.end))

        # Here we calculate the intersection between the two sets, and
        # if it's larger than our excessive size, we know that they're
        # overlapped
        ix = cas.intersection(cbs)
        if len(ix) >= excessive:
            bad += 1
            qc_features.append(gen_qc_feature(
                min(ix),
                max(ix),
                "Excessive Overlap")
            )
            results.append((gene_a, gene_b, min(ix), max(ix)))

    # Good isn't accurate here. It's a triangle number and just ugly, but we
    # don't care enough to fix it.
    return len(list(genes(record.features))), bad, results, qc_features


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
    return good, bad, results, []


def evaluate_and_report(annotations, genome, gff3):
    """
    Generate our HTML evaluation of the genome
    """
    # Get features from GFF file
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    # Get the first GFF3 record
    record = list(GFF.parse(annotations, base_dict=seq_dict))[0]

    gff3_qc_record = copy.deepcopy(record)
    gff3_qc_record.features = []
    gff3_qc_features = []
    upstream_min = 5
    upstream_max = 15

    log.info("Locating missing RBSs")
    mb_good, mb_bad, mb_results, mb_annotations = missing_rbs(
        record,
        lookahead_min=upstream_min,
        lookahead_max=upstream_max
    )
    gff3_qc_features += mb_annotations
    #log.info('%s %s %s', mb_good, mb_bad, mb_results)

    log.info("Locating excessive gaps")
    eg_good, eg_bad, eg_results, eg_annotations = excessive_gap(record, excess=3 * upstream_max)
    gff3_qc_features += eg_annotations

    log.info("Locating excessive overlaps")
    eo_good, eo_bad, eo_results, eo_annotations = excessive_overlap(record, excessive=15)
    gff3_qc_features += eo_annotations

    log.info("Locating morons")
    mo_good, mo_bad, mo_results, mo_annotations = find_morons(record)
    gff3_qc_features += mo_annotations

    score_good = float(sum((mb_good, eg_good, eo_good)))
    score_bad = float(sum((mb_bad, eg_bad, eo_bad)))

    score = int(100 * score_good / (score_bad + score_good))

    # This is data that will go into our HTML template
    kwargs = {
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

    with open(gff3, 'w') as handle:
        gff3_qc_record.features = gff3_qc_features
        gff3_qc_record.annotations = {}
        GFF.write([gff3_qc_record], handle)

    return REPORT_TEMPLATE.render(**kwargs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('annotations', type=file, help='Parent GFF3 annotations')
    parser.add_argument('genome', type=file, help='Genome Sequence')
    parser.add_argument('--gff3', type=str, help='GFF3 Annotations', default='qc_annotations.gff3')
    args = parser.parse_args()

    print evaluate_and_report(**vars(args))
