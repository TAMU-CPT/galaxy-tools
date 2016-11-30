#!/usr/bin/env python3.4
''' Parse Blast output in XML with Biopython and converts to SAM (v1).
Tested with Biopython 1.64 and BLASTN 2.2.30+ command

    blastn -task blastn -subject ref.fasta -query reads.fasta -outfmt 5 \
            -out outblast.xml -word_size 7 -qcov_hsp_perc 0.3

There are m times n records in blast xml output file, where m is the number of
sequences in the database (references) and n the number of queries (reads).

record -> alignment -> hsp

Record section
--------------

{'alignments': [<Bio.Blast.Record.Alignment object at 0x7feebecda9b0>],
 'application': 'BLASTN',
 'blast_cutoff': (None, None),
 'database': '',
 'database_length': 0,
 'database_letters': None,
 'database_name': [],
 'database_sequences': 0,
 'date': '',
 'descriptions': [<Bio.Blast.Record.Description object at 0x7feebecdaa58>],
 'dropoff_1st_pass': (None, None),
 'effective_database_length': None,
 'effective_hsp_length': 4,
 'effective_query_length': None,
 'effective_search_space': 621.0,
 'effective_search_space_used': None,
 'expect': '10',
 'filter': 'L;m;',
 'frameshift': (None, None),
 'gap_penalties': (5, 2),
 'gap_trigger': (None, None),
 'gap_x_dropoff': (None, None),
 'gap_x_dropoff_final': (None, None),
 'gapped': 0,
 'hsps_gapped': None,
 'hsps_no_gap': None,
 'hsps_prelim_gapped': None,
 'hsps_prelim_gapped_attemped': None,
 'ka_params': (0.625, 0.41, 0.78),
 'ka_params_gap': (None, None, None),
 'matrix': '',
 'multiple_alignment': None,
 'num_good_extends': None,
 'num_hits': None,
 'num_letters_in_database': 0,
 'num_seqs_better_e': None,
 'num_sequences': None,
 'num_sequences_in_database': 0,
 'posted_date': [],
 'query': '1X',
 'query_id': 'Query_1',
 'query_length': 13,
 'query_letters': 13,
 'reference': 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. '
              'Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and '
              'David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new '
              'generation of protein database search programs", Nucleic '
              'Acids Res. 25:3389-3402.',
 'sc_match': 2,
 'sc_mismatch': -3,
 'threshold': None,
 'version': '2.2.30+',
 'window_size': None}

Alignment section
-----------------

{'accession': 'Subject_1',
 'hit_def': 'mock_ref',
 'hit_id': 'Subject_1',
 'hsps': [<Bio.Blast.Record.HSP object at 0x7f4034696ba8>],
 'length': 73,
 'title': 'Subject_1 mock_ref'}

HSP section
-----------

hsps is a list, one hsp is
{'align_length': 13,
 'bits': 19.32,
 'expect': 0.000948843,
 'frame': (1, 1),
 'gaps': 0,
 'identities': 12,
 'match': '||| |||||||||',
 'num_alignments': None,
 'positives': 12,
 'query': 'GACAGATTACAGT',
 'query_end': 13,
 'query_start': 1,
 'sbjct': 'GACTGATTACAGT',
 'sbjct_end': 64,
 'sbjct_start': 52,
 'score': 20.0,
 'strand': (None, None)}



SAM alignment mandatory fields
1 QNAME String [!-?A-~]{1, 255} Query template NAME
2 FLAG Int [0, 2^16 - 1] bitwise FLAG
3 RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
4 POS Int [0, 2^31 - 1] 1-based leftmost mapping POSition
5 MAPQ Int [0, 255] MAPping Quality
6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next read
8 PNEXT Int [0,2^31 - 1] Position of the mate/next read
9 TLEN Int [-2^31 + 1, 2^31 - 1] observed Template LENgth
10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity + 33
'''

import sys
import copy
from itertools import tee
from math import log10
from Bio.Blast import NCBIXML

NOMAPQ = False
def_qual = 'I'

# Usage
try:
    filein = sys.argv[1]
except KeyError:
    sys.exit("Usage: blast2sam2.py <in+blastn>\n")


sam_line = ['', 0, None, 0, 0, None, '*', 0, 0, '*', '*']


def cigar(subject, query, queryStart, queryEnd, querySize):
    '''Build CIGAR representation from an HSP

    GTCCATGCAATTTTAAGACTTGAACCCCCTTGACTGATTACAGTCAGT   original sequence: 48 bp

          22 matches    +    14 gaps      +     26 matches  =  querySize: 62 bp
    GTCCATGCAATTTTAAGACTTG--------------AACCCCCTTGACTGATTACAGTCAGT  query
    ||||||||||||||||||||||              ||||||||||||||||||||||||||  midline
    GTCCATGCAATTTTAAGACTTGAACCTGTGATCTGAAACCCCCTTGACTGATTACAGTCAGT  subject
    '''

    # To store CIGAR representation
    cigar_str = []
    # Head clipping
    if queryStart > 1:
        cigar_str.append('%dH' % (queryStart - 1))

    # Evaluate alignment position by position (always begin with a match)
    curType = "="
    prevType = "="
    count = 0
    cigarsum = 0

    assert len(query) == querySize
    length = len(query.replace('-', ''))

    # this loops over the alignment positions
    for i in range(querySize):
        # Current position type (deletion, insertion, match or mismatch)
        if query[i] == '-':
            curType = 'D'
        elif subject[i] == '-':
            curType = 'I'
        elif query[i] == subject[i]:
            curType = '='
        else:
            curType = 'X'

        if curType == prevType:
            # Enlarge current segment
            prevType = curType
            count += 1
        else:
            # Write current segment and start a new one
            cigar_str.append('%d%s' % (count, prevType))
            if prevType in ['I', '=', 'X']:
                cigarsum += count
            prevType = curType
            count = 1

    # Write last group
    cigar_str.append('%d%s' % (count, curType))
    if curType in ['I', '=', 'X']:
        cigarsum += count

    # Tail clipping
    if queryEnd < length:
        # print(query, querySize, queryEnd, file=sys.stderr)
        cigar_str.append('%dH' % (length - queryEnd))

    assert cigarsum == length, '%s: cigar:%s\tcigarsum=%d,length=%d' % \
        (query, ''.join(cigar_str), cigarsum, length)
    # Join segments into a string
    return ''.join(cigar_str)


# use itertools.tee() because we need the list twice
blast_records, blast_records_backup = tee(NCBIXML.parse(open(filein)))

# read once to parse general info
version = None
references = {}
for record in blast_records:
    if not version:
        version = record.version
        application = record.application
    if record.alignments == []:
        continue
    for alignment in record.alignments:
        # references[alignment.hit_def] = alignment.length
        references[record.query] = alignment.length

# print header
print('@HD\tVN:1.0\tSO:unsorted')
for k, v in references.items():
    print('@SQ\tSN:%s\tLN:%d' % (k, v))
print('@PG\tID:%s\tVN:%s\tCL:%s' % (application, version, ' '.join(sys.argv)))

counter = {}
i = 0
for record in blast_records_backup:
    for alignment in record.alignments:
        TC = len(alignment.hsps)  # SAM TC flag: segments in template
        for hsp in alignment.hsps:
            to_print = copy.copy(sam_line)

            idx = 0
            k = '|||'.join((record.query, alignment.hit_id))
            if k in counter:
                counter[k] += 1
                idx = counter[k]
            else:
                counter[k] = 0

            to_print[0] = alignment.hit_id + ':%s' % idx
            i += 1
            to_print[2] = record.query
            # to_print[0] = record.query
            # to_print[2] = alignment.hit_def
            to_print[3] = min(hsp.sbjct_start, hsp.sbjct_end)

            try:
                mapq = int(-log10(hsp.expect))
            except ValueError:
                mapq = 255

            if mapq > 254:
                to_print[4] = 254
            elif mapq < 0:
                to_print[4] = 0
            else:
                to_print[4] = mapq
            if NOMAPQ:
                to_print[4] = 255

            if hsp.frame == (1, -1):
                to_print[1] |= 16
            elif hsp.frame == (1, 1):
                pass
            else:
                # Handle in sam_rebase script
                pass

            to_print[5] = cigar(hsp.sbjct, hsp.query, hsp.query_start,
                                hsp.query_end, hsp.align_length)

            to_print[9] = hsp.query.replace('-', '')
            to_print[10] = def_qual * len(to_print[9])
            to_print = [str(t) for t in to_print]
            print('\t'.join(to_print))
