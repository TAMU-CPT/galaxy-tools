#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import sys
import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def get_id(feature=None, parent_prefix=None):
    result = ""
    if parent_prefix is not None:
        result += parent_prefix + "|"
    if "locus_tag" in feature.qualifiers:
        result += feature.qualifiers["locus_tag"][0]
    elif "gene" in feature.qualifiers:
        result += feature.qualifiers["gene"][0]
    elif "product" in feature.qualifiers:
        result += feature.qualifiers["product"][0]
    else:
        result += "%s_%s_%s" % (
            feature.location.start,
            feature.location.end,
            feature.location.strand,
        )
    return result


def ensure_location_in_bounds(start=0, end=0, parent_length=0):
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


def extract_features(
    genbank_file=None,
    tag="CDS",
    translate=False,
    n_bases_upstream=0,
    n_bases_downstream=0,
    strip_stops=False,
    translation_table_id=11,
    informative=False,
):

    codon_table = CodonTable.unambiguous_dna_by_id[11]
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type in tag:
              finSeq = ""
              for locInd in range(0, len(feature.location.parts)):
                thisLoc = feature.location.parts[locInd]
                # Find new feature boundaries
                start = int(thisLoc.start)
                end = int(thisLoc.end)
                strand = thisLoc.strand
                if strand < 0:
                  temp = end
                  end = min(start, end)
                  start = max(start, temp)
                
                if locInd == 1:
                  start += 0

                if n_bases_downstream != 0:
                    # If we want extra on the end we cannot listen to
                    # stop_stripping requests
                    if strand >= 0 and locInd == len(feature.location.parts) - 1:
                        end += n_bases_downstream
                    elif strand < 0 and locInd == 0:
                        start -= n_bases_downstream

                # n_bases_upstream
                if strand >= 0 and locInd == 0:
                    start -= n_bases_upstream
                elif strand < 0 and locInd == len(feature.location.parts) - 1:
                    end += n_bases_upstream

                __seqs = []
                # Upstream addition
                if n_bases_upstream > 0:
                    __seqs.append(
                        SeqFeature(
                            FeatureLocation(
                                start, int(thisLoc.start), strand=strand
                            ),
                            type="domain",
                        )
                    )

                __seqs.append(SeqFeature(
                                FeatureLocation(
                                start, end, strand=strand
                                ),
                                type="domain",
                             ))
                # Downstream addition
                if n_bases_downstream > 0:
                    __seqs.append(
                        SeqFeature(
                            FeatureLocation(
                                int(thisLoc.end), end, strand=strand
                            ),
                            type="domain",
                        )
                    )

                rangeS = __seqs[0].location.start
                rangeE = __seqs[0].location.end
                for s in __seqs:
                    if rangeS > s.location.start:
                        rangeS = s.location.start
                    if rangeE < s.location.end:
                        rangeE = s.location.end

                if "codon_start" in feature.qualifiers:
                    if strand > 0:
                        rangeS += int(feature.qualifiers["codon_start"][0]) - 1
                    else:
                        rangeE -= int(feature.qualifiers["codon_start"][0]) - 1

                if translate:
                    try:
                        if strand > 0:
                            retSeq = (record.seq[rangeS:rangeE]).translate(
                                table=translation_table_id, cds=True
                            )
                        else:
                            retSeq = (
                                (record.seq[rangeS:rangeE])
                                .reverse_complement()
                                .translate(table=translation_table_id, cds=True)
                            )
                    except Exception as bdct:
                        if strand > 0:
                            retSeq = (record.seq[rangeS:rangeE]).translate(
                                table=translation_table_id, cds=False
                            )
                        else:
                            retSeq = (
                                (record.seq[rangeS:rangeE])
                                .reverse_complement()
                                .translate(table=translation_table_id, cds=False)
                            )
                            log.warn("ERROR %s %s", get_id(feature), bdct)
                    # extracted_seqs = []
                    # for x in __seqs:
                    #    try:
                    #        y = x.extract(record.seq).translate(table=translation_table_id, cds=True)
                    #        extracted_seqs.append(y)
                    #    except Exception as bdct:
                    #        y = x.extract(record.seq).translate(table=translation_table_id, cds=False)
                    #        extracted_seqs.append(y)
                    #        log.warn("ERROR %s %s", get_id(x), bdct)
                else:
                    if strand > 0:
                            retSeq = (record.seq[rangeS:rangeE])
                    else:
                            retSeq = (
                                (record.seq[rangeS:rangeE])
                                .reverse_complement())
                    # extracted_seqs = [x.extract(record.seq) for x in __seqs]

                if informative:
                  if locInd == 0:
                    defline = " %s [start=%s,end=%s]" % (
                        ",".join(feature.qualifiers.get("product", [])),
                        start + 1,
                        end,
                    )
                  else:
                    defline += ", joined with [start=%s,end=%s]" % (start + 1, end)
                   
                else:
                  if locInd == 0:
                    defline = " [start=%s,end=%s]" % (start + 1, end)
                  else:
                    defline += ", joined with [start=%s,end=%s]" % (start + 1, end)

                extracted_seq = str(retSeq)
                
                if strip_stops:
                    if extracted_seq[-3:] in codon_table.stop_codons:
                        extracted_seq = extracted_seq[:-3]
                    elif extracted_seq[-1:] in "*":
                        extracted_seq = extracted_seq[:-1]
                finSeq += extracted_seq

              yield [
                    SeqRecord(
                        Seq(finSeq.strip()),
                        id=get_id(feature),
                        description=defline,
                    )
              ]


if __name__ == "__main__":
    # Grab all of the filters from our plugin loader
    gbk_tags = [
        "all",
        "-10_signal",
        "-35_signal",
        "3'UTR",
        "5'UTR",
        "CAAT_signal",
        "CDS",
        "C_region",
        "D-loop",
        "D_segment",
        "GC_signal",
        "J_segment",
        "LTR",
        "N_region",
        "RBS",
        "STS",
        "S_region",
        "TATA_signal",
        "V_region",
        "V_segment",
        "assembly_gap",
        "attenuator",
        "enhancer",
        "exon",
        "gap",
        "gene",
        "iDNA",
        "intron",
        "mRNA",
        "mat_peptide",
        "misc_RNA",
        "misc_binding",
        "misc_difference",
        "misc_feature",
        "misc_recomb",
        "misc_signal",
        "misc_structure",
        "mobile_element",
        "modified_base",
        "ncRNA",
        "old_sequence",
        "operon",
        "oriT",
        "polyA_signal",
        "polyA_site",
        "precursor_RNA",
        "prim_transcript",
        "primer_bind",
        "promoter",
        "protein_bind",
        "rRNA",
        "rep_origin",
        "repeat_region",
        "sig_peptide",
        "source",
        "stem_loop",
        "tRNA",
        "terminator",
        "tmRNA",
        "transit_peptide",
        "unsure",
        "variation",
    ]

    parser = argparse.ArgumentParser(
        description="Export a subset of features from a Genbank file", epilog=""
    )
    parser.add_argument(
        "genbank_file", type=argparse.FileType("r"), help="Genbank file"
    )
    parser.add_argument(
        "tag", nargs="+", type=str, choices=gbk_tags, help="tags to export"
    )

    parser.add_argument("--translate", action="store_true", help="Translate sequence")
    parser.add_argument(
        "--translation_table_id", help="Translation table ID", default=11
    )
    parser.add_argument(
        "--n_bases_upstream",
        type=int,
        help="Add N bases upstream to exported features",
        default=0,
    )
    parser.add_argument(
        "--n_bases_downstream",
        type=int,
        help="Add N bases downstream to exported features",
        default=0,
    )
    parser.add_argument("--strip_stops", action="store_true", help="Remove stop codons")
    parser.add_argument(
        "--informative", action="store_true", help="More informative deflines"
    )

    args = vars(parser.parse_args())
    for seq in extract_features(**args):
        SeqIO.write(seq, sys.stdout, "fasta")
