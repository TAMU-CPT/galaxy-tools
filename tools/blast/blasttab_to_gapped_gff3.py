#!/usr/bin/env python
import re
import sys
import copy
import argparse
from BCBio import GFF
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name="blasttab2gff3")

__doc__ = """
Blast TSV files, when transformed to GFF3, do not normally show gaps in the
blast hits. This tool aims to fill that "gap".
"""


def blasttsv2gff3(
    blasttsv, min_gap=3, trim_start=False, trim_end=False, type="nucleotide_match"
):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    # http://www.sequenceontology.org/browser/release_2.4/term/SO:0000343
    match_type = {  # Currently we can only handle BLASTN, BLASTP
        "BLASTN": "nucleotide_match",
        "BLASTP": "protein_match",
    }.get(type, "match")

    columns = [
        "qseqid",  # 01 Query Seq-id (ID of your sequence)
        "sseqid",  # 02 Subject Seq-id (ID of the database hit)
        "pident",  # 03 Percentage of identical matches
        "length",  # 04 Alignment length
        "mismatch",  # 05 Number of mismatches
        "gapopen",  # 06 Number of gap openings
        "qstart",  # 07 Start of alignment in query
        "qend",  # 08 End of alignment in query
        "sstart",  # 09 Start of alignment in subject (database hit)
        "send",  # 10 End of alignment in subject (database hit)
        "evalue",  # 11 Expectation value (E-value)
        "bitscore",  # 12 Bit score
        "sallseqid",  # 13 All subject Seq-id(s), separated by a ';'
        "score",  # 14 Raw score
        "nident",  # 15 Number of identical matches
        "positive",  # 16 Number of positive-scoring matches
        "gaps",  # 17 Total number of gaps
        "ppos",  # 18 Percentage of positive-scoring matches
        "qframe",  # 19 Query frame
        "sframe",  # 20 Subject frame
        "qseq",  # 21 Aligned part of query sequence
        "sseq",  # 22 Aligned part of subject sequence
        "qlen",  # 23 Query sequence length
        "slen",  # 24 Subject sequence length
        "salltitles",  # 25 All subject title(s), separated by a '<>'
    ]

    for record_idx, record in enumerate(blasttsv):
        if record.startswith("#"):
            continue

        dc = {k: v for (k, v) in zip(columns, (x.strip() for x in record.split("\t")))}

        rec = SeqRecord(Seq("ACTG"), id=dc["qseqid"])

        feature_id = "blast.%s.%s.%s" % (record_idx, dc["qseqid"], dc["sseqid"])
        feature_id = re.sub("\|", "_", feature_id)
        feature_id = re.sub("[^A-Za-z0-9_.-]", "", feature_id)
        qualifiers = {
            "ID": feature_id,
            "Name": dc["salltitles"].split("<>")[0],
            "description": "Hit to {sstart}..{send} ({sframe}) of {x}".format(
                x=dc["salltitles"].split("<>")[0], **dc
            ),
            "source": "blast",
            "score": dc["evalue"],
            "accession": dc["sseqid"],
            "length": dc["qlen"],
            "hit_titles": dc["salltitles"].split("<>"),
            "Target": dc["qseqid"],
        }

        for key in dc.keys():
            if key in (
                "salltitles",
                "sallseqid",
                "score",
                "sseqid",
                "qseqid",
                "qseq",
                "sseq",
            ):
                continue
            qualifiers["blast_%s" % key] = dc[key]

        for (
            integer_numerical_key
        ) in "gapopen gaps length mismatch nident positive qend qframe qlen qstart score send sframe slen sstart".split(
            " "
        ):
            dc[integer_numerical_key] = int(dc[integer_numerical_key])

        for float_numerical_key in "bitscore evalue pident ppos".split(" "):
            dc[float_numerical_key] = float(dc[float_numerical_key])

        # This required a fair bit of sketching out/match to figure out
        # the first time.
        #
        # the match_start location must account for queries and
        # subjecst that start at locations other than 1
        parent_match_start = dc["qstart"]
        # The end is the start + hit.length because the match itself
        # may be longer than the parent feature, so we use the supplied
        # subject/hit length to calculate the real ending of the target
        # protein.
        parent_match_end = dc["qend"]

        # However, if the user requests that we trim the feature, then
        # we need to cut the ``match`` start to 0 to match the parent feature.
        # We'll also need to cut the end to match the query's end. It (maybe)
        # should be the feature end? But we don't have access to that data, so
        # We settle for this.

        # The ``match`` feature will hold one or more ``match_part``s
        top_feature = SeqFeature(
            FeatureLocation(
                min(parent_match_start, parent_match_end) - 1,
                max(parent_match_start, parent_match_end),
            ),
            type=match_type,
            strand=0,
            qualifiers=qualifiers,
        )
        top_feature.sub_features = []

        # Unlike the parent feature, ``match_part``s have sources.
        part_qualifiers = {"source": "blast"}
        for start, end, cigar in generate_parts(
            dc["qseq"], dc.get("mseq", None), dc["sseq"], ignore_under=min_gap
        ):

            part_qualifiers["Gap"] = cigar
            part_qualifiers["ID"] = dc["sseqid"]

            match_part_start = parent_match_start + start

            # We used to use hsp.align_length here, but that includes
            # gaps in the parent sequence
            #
            # Furthermore align_length will give calculation errors in weird places
            # So we just use (end-start) for simplicity
            match_part_end = match_part_start + (end - start)
            # print start, end, cigar, parent_match_start, parent_match_end, match_part_start, match_part_end, dc['qstart'], dc['qend'], dc['qframe'], dc['sstart'], dc['send'], dc['sframe'], dc['length'], dc['qseq'].count('-')

            top_feature.sub_features.append(
                SeqFeature(
                    FeatureLocation(
                        min(match_part_start, match_part_end) - 1,
                        max(match_part_start, match_part_end) - 1,
                    ),
                    type="match_part",
                    strand=0,
                    qualifiers=copy.deepcopy(part_qualifiers),
                )
            )

        top_feature.sub_features = sorted(
            top_feature.sub_features, key=lambda x: int(x.location.start)
        )
        rec.features = [top_feature]
        yield rec


def __remove_query_gaps(query, match, subject):
    """remove positions in all three based on gaps in query

    In order to simplify math and calculations...we remove all of the gaps
    based on gap locations in the query sequence::

        Q:ACTG-ACTGACTG
        S:ACTGAAC---CTG

    will become::

        Q:ACTGACTGACTG
        S:ACTGAC---CTG

    which greatly simplifies the process of identifying the correct location
    for a match_part
    """
    prev = 0
    fq = ""
    fm = ""
    fs = ""
    for position in re.finditer("-", query):
        fq += query[prev : position.start()]
        if match is not None:
            fm += match[prev : position.start()]
        fs += subject[prev : position.start()]
        prev = position.start() + 1
    fq += query[prev:]
    if match is not None:
        fm += match[prev:]
    fs += subject[prev:]

    return (fq, None if match is None else fm, fs)


def _none_safe_match_iterator(match):
    if match is None:
        while True:
            yield None
    else:
        for x in match:
            yield x


def generate_parts(query, match, subject, ignore_under=3):
    region_q = []
    region_m = []
    region_s = []

    (query, match, subject) = __remove_query_gaps(query, match, subject)

    region_start = -1
    region_end = -1
    mismatch_count = 0
    for i, (q, m, s) in enumerate(
        zip(query, _none_safe_match_iterator(match), subject)
    ):

        # If we have a match
        if (m is None and q == s) or (m is not None and (m != " " or m == "+")):
            if region_start == -1:
                region_start = i
                # It's a new region, we need to reset or it's pre-seeded with
                # spaces
                region_q = []
                region_m = []
                region_s = []
            region_end = i
            mismatch_count = 0
        else:
            mismatch_count += 1

        region_q.append(q)
        region_m.append(m)
        region_s.append(s)

        if mismatch_count >= ignore_under and region_start != -1 and region_end != -1:
            region_q = region_q[0:-ignore_under]
            region_m = region_m[0:-ignore_under]
            region_s = region_s[0:-ignore_under]
            yield region_start, region_end + 1, cigar_from_string(
                region_q, region_m, region_s, strict_m=True
            )
            region_q = []
            region_m = []
            region_s = []

            region_start = -1
            region_end = -1
            mismatch_count = 0

    yield region_start, region_end + 1, cigar_from_string(
        region_q, region_m, region_s, strict_m=True
    )


def _qms_to_matches(query, match, subject, strict_m=True):
    matchline = []

    for (q, m, s) in zip(query, match, subject):
        ret = ""

        if m != " " or m == "+":
            ret = "="
        elif m == " ":
            if q == "-":
                ret = "D"
            elif s == "-":
                ret = "I"
            else:
                ret = "X"
        else:
            log.warn("Bad data: \n\t%s\n\t%s\n\t%s\n" % (query, match, subject))

        if strict_m:
            if ret == "=" or ret == "X":
                ret = "M"

        matchline.append(ret)
    return matchline


def _matchline_to_cigar(matchline):
    cigar_line = []
    last_char = matchline[0]
    count = 0
    for char in matchline:
        if char == last_char:
            count += 1
        else:
            cigar_line.append("%s%s" % (last_char, count))
            count = 1
        last_char = char
    cigar_line.append("%s%s" % (last_char, count))
    return " ".join(cigar_line)


def cigar_from_string(query, match, subject, strict_m=True):
    matchline = _qms_to_matches(query, match, subject, strict_m=strict_m)
    if len(matchline) > 0:
        return _matchline_to_cigar(matchline)
    else:
        return ""


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Blast TSV to gapped GFF3")
    parser.add_argument(
        "blasttsv", type=argparse.FileType("r"), help="Blast TSV Output"
    )
    parser.add_argument(
        "--min_gap",
        type=int,
        help="Maximum gap size before generating a new match_part",
        default=3,
    )
    parser.add_argument(
        "--trim_start",
        action="store_true",
        help="Trim blast hits to be only as long as the parent feature",
    )
    parser.add_argument(
        "--trim_end", action="store_true", help="Cut blast results off at end of gene"
    )
    parser.add_argument("--type", choices=["BLASTN"])
    args = parser.parse_args()

    for rec in blasttsv2gff3(**vars(args)):
        GFF.write([rec], sys.stdout)
