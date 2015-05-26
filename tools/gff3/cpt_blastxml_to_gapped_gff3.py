#!/usr/bin/perl
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import sys
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name='blastxml2gff3')
from BCBio import GFF


__doc__ = """
BlastXML files, when transformed to GFF3, do not normally show gaps in the
blast hits. This tool aims to fill that "gap".

For an input BlastXML file, this tool will produce a GFF3 file containing all
of the relevant information: start, stop, score, and the GAP as a CIGAR string.

`One CIGAR Spec <http://www.sequenceontology.org/gff3.shtml>`__ does not list X
and = for mismatches and matches, respsectively. However,
`Other CIGAR Specs <http://samtools.github.io/hts-specs/SAMv1.pdf>`__ do allow
for those characters. Given that I cannot anticipate which characters your
downstream analysis methods will support, this tool provides the option to use
just M (strict_m enabled), or to use both, as is the default.
"""

def blastxml2gff3(blastxml=None, strict_m=False, **kwd):
    from Bio.Blast import NCBIXML
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    blast_records = NCBIXML.parse(blastxml)
    records = []
    for record in blast_records:
        #'alignments', 'application', 'blast_cutoff', 'database',
        #'database_length', 'database_letters', 'database_name',
        #'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass',
        #'effective_database_length', 'effective_hsp_length',
        #'effective_query_length', 'effective_search_space',
        #'effective_search_space_used', 'expect', 'filter', 'frameshift',
        #'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final',
        #'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped',
        #'hsps_prelim_gapped_attemped', 'ka_params', 'ka_params_gap', 'matrix',
        #'multiple_alignment', 'num_good_extends', 'num_hits',
        #'num_letters_in_database', 'num_seqs_better_e', 'num_sequences',
        #'num_sequences_in_database', 'posted_date', 'query', 'query_id',
        #'query_length', 'query_letters', 'reference', 'sc_match',
        #'sc_mismatch', 'threshold', 'version', 'window_size'

        rec = SeqRecord(Seq("ACTG"), id=record.query)
        for hit in record.alignments:
            #'accession', 'hit_def', 'hit_id', 'hsps', 'length', 'title'
            for hsp in hit.hsps:
                #'align_length', 'bits', 'expect', 'frame', 'gaps',
                #'identities', 'match', 'num_alignments', 'positives', 'query',
                #'query_end', 'query_start', 'sbjct', 'sbjct_end',
                #'sbjct_start', 'score', 'strand'
                qualifiers = {"source": "blast",
                              "score": hsp.expect,
                              "accession": hit.accession,
                              "hit_id": hit.hit_id,
                              "length": hit.length,
                              "title": hit.title,
                              "Gap": cigar_from_string(hsp.query, hsp.match, hsp.sbjct, strict_m)}
                top_feature = SeqFeature(FeatureLocation(hsp.sbjct_start,
                                                         hsp.sbjct_end),
                                         type="Match", strand=0,
                                         qualifiers=qualifiers)
                rec.features.append(top_feature)
        records.append(rec)
    return records


def cigar_from_string(query, match, subject, strict_m=True):
    last = None
    last_count = 0
    cigar_line = ''

    for (q, m, s) in zip(query, match, subject):
        ret = ''

        if m != ' ' or m == '+':
            ret = '='
        elif m == ' ':
            if q == '-':
                ret = 'D'
            elif s == '-':
                ret = 'I'
            else:
                ret = 'X'
        else:
            log.warn("Bad data: \n\t%s\n\t%s\n\t%s\n" % (query, match, subject))

        if strict_m:
            if ret == '=' or ret == 'X':
                ret = 'M'

        last_count += 1
        if ret != last:
            if last is not None:
                cigar_line += "%s%s " % (last, last_count)
            last_count = 0
            last = ret
    if last is not None:
        cigar_line += "%s%s " % (last, last_count)

    return cigar_line


def passthrough(cb):
    opts = GGO(
        options=[
            ['blastxml', 'Input Blast XML Data',
             {'required': True, 'validate': 'File/Input'}],
            [ 'strict_m', 'One CIGAR spec specificies that Matches AND Mismatches are both represented as M. Another spec allows for = and X to disambiguate. Depending on your downstream pipeline, choose to use =/X or not.',
             { 'validate' : 'Flag' }],
        ],
        outputs=[
            [
                'gff',
                'GFF Output File',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'gff',
                    'data_format': 'text/plain',
                    'default_format': 'TXT',
                }
            ],
        ],
        defaults={
            'appid': 'edu.tamu.cpt.blast.gapped_gff3',
            'appname': 'BlastXML to GFF3',
            'appvers': '1.0',
            'appdesc': 'converts blast XML data to gapped alignment',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    return (opts, options, cb(**options))


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    (opts, options, result) = passthrough(blastxml2gff3)

    with open(options['gff'], "w") as out_handle:
        GFF.write(result, out_handle)
