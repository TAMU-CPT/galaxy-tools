#!/usr/bin/env python
import argparse
import copy
import logging
import re
import sys
from BCBio import GFF

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name="blast2gff3")

__doc__ = """
BlastXML files, when transformed to GFF3, do not normally show gaps in the
blast hits. This tool aims to fill that "gap".
"""


def blast2gff3(blast, blastxml=False, blasttab=False, include_seq=False):
    # Call correct function based on xml or tabular file input, raise error if neither or both are provided
    if blastxml and blasttab:
        raise Exception("Cannot provide both blast xml and tabular file")

    if blastxml:
        recs = blastxml2gff3(blast, include_seq)
        return recs
    elif blasttab:
        raise Exception("Tabular not supported yet, please use xml")
        # yield blasttab2gff3(blast, include_seq)
    else:
        raise Exception("Must provide either blast xml or tabular file")


def blastxml2gff3(blastxml, include_seq=False):
    from Bio.Blast import NCBIXML
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    blast_records = NCBIXML.parse(blastxml)
    for idx_record, record in enumerate(blast_records):
        # http://www.sequenceontology.org/browser/release_2.4/term/SO:0000343
        match_type = {  # Currently we can only handle BLASTN, BLASTP
            "BLASTN": "nucleotide_match",
            "BLASTP": "protein_match",
        }.get(record.application, "match")

        recid = record.query
        if " " in recid:
            recid = recid[0 : recid.index(" ")]

        rec = SeqRecord(Seq("ACTG"), id=recid)
        for idx_hit, hit in enumerate(record.alignments):
            # gotta check all hsps in a hit to see boundaries
            parent_match_start = 0
            parent_match_end = 0
            hit_qualifiers = {
                "ID": "b2g.%s.%s.%s" % (idx_record, idx_hit, "0"),
                "source": "blast",
                "accession": hit.accession,
                "hit_id": hit.hit_id,
                "length": hit.length,
                "hit_titles": hit.title.split(" >"),
                "hsp_count": len(hit.hsps),
            }
            sub_features = []
            for idx_hsp, hsp in enumerate(hit.hsps):
                hsp_qualifiers = {
                    "ID": "b2g.%s.%s.%s" % (idx_record, idx_hit, idx_hsp),
                    "source": "blast",
                    "score": hsp.expect,
                    "accession": hit.accession,
                    "hit_id": hit.hit_id,
                    "length": hit.length,
                    "hit_titles": hit.title.split(" >"),
                }
                if include_seq:
                    hsp_qualifiers.update(
                        {
                            "blast_qseq": hsp.query,
                            "blast_sseq": hsp.sbjct,
                            "blast_mseq": hsp.match,
                        }
                    )

                for prop in (
                    "score",
                    "bits",
                    "identities",
                    "positives",
                    "gaps",
                    "align_length",
                    "strand",
                    "frame",
                    "query_start",
                    "query_end",
                    "sbjct_start",
                    "sbjct_end",
                ):
                    hsp_qualifiers["blast_" + prop] = getattr(hsp, prop, None)

                desc = hit.title.split(" >")[0]
                hsp_qualifiers["description"] = desc[desc.index(" ") :]

                # check if parent boundary needs to increase
                if hsp.query_start < parent_match_start:
                    parent_match_start = hsp.query_start
                if hsp.query_end > parent_match_end:
                    parent_match_end = hsp.query_end + 1

                # Build out the match_part features for each HSP
                for idx_part, (start, end, cigar) in enumerate(
                    generate_parts(hsp.query, hsp.match, hsp.sbjct, ignore_under=10)
                ):
                    hsp_qualifiers["Gap"] = cigar
                    hsp_qualifiers["ID"] = hit_qualifiers["ID"] + (".%s" % idx_part)

                    match_part_start = hsp.query_start

                    # We used to use hsp.align_length here, but that includes
                    # gaps in the parent sequence
                    #
                    # Furthermore align_length will give calculation errors in weird places
                    # So we just use (end-start) for simplicity
                    match_part_end = match_part_start + (end - start)

                    sub_features.append(
                        SeqFeature(
                            FeatureLocation(match_part_start, match_part_end),
                            type="match_part",
                            strand=0,
                            qualifiers=copy.deepcopy(hsp_qualifiers),
                        )
                    )

            # Build the top level seq feature for the hit
            top_feature = SeqFeature(
                FeatureLocation(parent_match_start, parent_match_end),
                type=match_type,
                strand=0,
                qualifiers=hit_qualifiers,
            )
            # add the generated subfeature hsp match_parts to the hit feature
            top_feature.sub_features = copy.deepcopy(sub_features)
            # Add the hit feature to the record
            rec.features.append(top_feature)
        rec.annotations = {}
        yield rec


def blasttab2gff3(blasttab, include_seq=False):
    yield none


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
        fm += match[prev : position.start()]
        fs += subject[prev : position.start()]
        prev = position.start() + 1
    fq += query[prev:]
    fm += match[prev:]
    fs += subject[prev:]

    return (fq, fm, fs)


def generate_parts(query, match, subject, ignore_under=3):
    region_q = []
    region_m = []
    region_s = []

    (query, match, subject) = __remove_query_gaps(query, match, subject)

    region_start = -1
    region_end = -1
    mismatch_count = 0
    for i, (q, m, s) in enumerate(zip(query, match, subject)):

        # If we have a match
        if m != " " or m == "+":
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
    parser = argparse.ArgumentParser(
        description="Convert Blast output to gapped GFF3, must provide XML or Tabular output",
        epilog="",
    )
    parser.add_argument(
        "blast",
        type=argparse.FileType("r"),
        help="Blast XML or 25 Column Tabular Output file",
    )
    parser.add_argument(
        "--blastxml", action="store_true", help="Process file as Blast XML Output"
    )
    parser.add_argument(
        "--blasttab",
        action="store_true",
        help="Process file as Blast 25 Column  Tabular Output",
    )
    parser.add_argument("--include_seq", action="store_true", help="Include sequence")
    args = parser.parse_args()

    for rec in blast2gff3(**vars(args)):
        if len(rec.features):
            GFF.write([rec], sys.stdout)
