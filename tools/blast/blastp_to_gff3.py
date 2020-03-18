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
        #yield blasttab2gff3(blast, include_seq)
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
                    "ID": "b2g.%s.%s.%s" % (idx_record, idx_hit, '0'),
                    "source": "blast",
                    "accession": hit.accession,
                    "hit_id": hit.hit_id,
                    "length": hit.length,
                    "hit_titles": hit.title.split(" >"),
                    "hsp_count": len(hit.hsps),
                    }
            sub_features = []
            for idx_hsp, hsp in enumerate(hit.hsps):
                if idx_hsp == 0:
                    # -1 and +1 for start/end to convert 0 index of python to 1 index of people
                    parent_match_start = hsp.query_start - 1
                    parent_match_end = hsp.query_end + 1
                # generate qualifiers to be added to gff3 feature
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

                #check if parent boundary needs to increase to envelope hsp
                if hsp.query_start < parent_match_start:
                    parent_match_start = hsp.query_start - 1
                if hsp.query_end > parent_match_end:    
                    parent_match_end = hsp.query_end + 1

                # add hsp to the gff3 feature as a "match_part"
                sub_features.append(
                    SeqFeature(
                        FeatureLocation(hsp.query_start, hsp.query_end),
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


def blasttsv2gff3(
    blasttsv, include_seq=False 
):
    yield None

def blasttab2gff3(blasttab, include_seq=False):
    yield None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert Blast output to gapped GFF3, must provide XML or Tabular output", epilog=""
    )
    parser.add_argument(
            "blast", type=argparse.FileType("r"), help="Blast XML or 25 Column Tabular Output file"
    )
    parser.add_argument(
        "--blastxml", action="store_true", help="Process file as Blast XML Output"
    )
    parser.add_argument(
        "--blasttab", action="store_true", help="Process file as Blast 25 Column  Tabular Output"
    )
    parser.add_argument("--include_seq", action="store_true", help="Include sequence")
    args = parser.parse_args()

    for rec in blast2gff3(**vars(args)):
        if len(rec.features):
            GFF.write([rec], sys.stdout)
