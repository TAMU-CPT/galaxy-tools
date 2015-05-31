#!/usr/bin/env python
import re
import sys
import argparse
from Bio import SeqIO

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()
LIPOBOX = re.compile('[ILMFTV][^REKD][GAS]C')

"""
1.  Take all the NONEs and sort out all those with genome size <90 kb.
2.  Call this the SMALLNONEs.  Everything else is PROBLEMATIC (because
    lysis cassettes are not as common in large genomes)
3.  Then sort out all the genomes with the words holin and endolysin anywhere in the record.
    a.   If I am right, holin and endolysin will be within 3 CDSs of each
         other.  If not, then put into PROBLEMATIC
    b.   Then look for all possible Lipoboxes within x= 2kb of the endolysin
         CDS. Let x and 'endolysin' be default but assignable.
    c.   Then manually inspect each lipobox possibility ....
"""

"""
3.  Then sort out all the genomes with the words holin and endolysin anywhere in the record.
    a.   If I am right, holin and endolysin will be within 3 CDSs of each
         other.  If not, then put into PROBLEMATIC
    b.   Then look for all possible Lipoboxes within x= 2kb of the endolysin
         CDS. Let x and 'endolysin' be default but assignable.
    c.   Then manually inspect each lipobox possibility ....
"""


__BASE_ENDOLYSIN_KEYWORDS = (
    "AMI-2",
    "AMI-3",
    "AMI-5",
    "AMI02-C",
    "Amidase02_C ",
    "Amidase_2  ",
    "Amidase_3 ",
    "Amidase_5 ",
    "CHAP",
    "CPL7",
    "ChW",
    "ChW",
    "Chitinase",
    "Cpl-7",
    "DUF3597",
    "GH108 ",
    "GH19",
    "GH25",
    "Glyco_hydro_108",
    "Glyco_hydro_19 ",
    "Glyco_hydro_25",
    "GyH",
    "LGFP",
    "LysM",
    "N-acetylglucosaminidase",
    "N-acetylmuramidase",
    "N-acetylmuramidase",
    "N-acetylmuraminidase",
    "NLPC_P60",
    "NLPC-P60",
    "NlpD",
    "PET-15-3",
    "PET-15-4",
    "PET-C39-2",
    "PET-M23",
    "PET-U40",
    "PG-1",
    "PG-3",
    "PG_binding_1",
    "PG_binding_3",
    "Peptidase_C39_2",
    "Peptidase_M15_3",
    "Peptidase_M15_4",
    "Peptidase_M23",
    "Peptidase_U40",
    "SH3-3",
    "SH3-5",
    "SH3-r",
    "SH3-related FOG",
    "SH3_3",
    "SH3_5",
    "SH3b ",
    "SLAP",
    "SLT",
    "SPOR",
    "Transglycosylase Glucosaminidase",
    "VanY",
    "YkuD",
    "amidase",
    "endolysin",
    "endopeptidase",
    "lysin",
    "lysozyme",
    "muramidase",
    "muraminidase",
    "peptidase",
)
ENDOLYSIN_KEYWORDS = []
for x in __BASE_ENDOLYSIN_KEYWORDS:
    ENDOLYSIN_KEYWORDS.append(x.lower())
    ENDOLYSIN_KEYWORDS.append(x.lower().replace('-', '_'))
    ENDOLYSIN_KEYWORDS.append(x.lower().replace('_', '-'))

HOLIN_KEYWORDS = (
    'holin'
)


def contains_keyword(feats, KWD_LIST):
    for feat in feats:
        for qual_key in feat.qualifiers.keys():
            for qual_val in feat.qualifiers[qual_key]:
                for kwd in KWD_LIST:
                    if kwd in qual_val:
                        yield feat


def feat_contains_keyword(feat, KWD_LIST):
    for qual_key in feat.qualifiers.keys():
        for qual_val in feat.qualifiers[qual_key]:
            for kwd in KWD_LIST:
                if kwd in qual_val:
                    return True
    return False


def features_in_range(features, center_feature, dist=2000):
    for feature in features:
        if feature.start < center_feature.start:
            if abs(center_feature.start - feature.end) < dist:
                yield feature
        else:
            if abs(center_feature.end - feature.start) < dist:
                yield feature


def tag(feat, tag_name, tag_value):
    if tag_name in feat:
        feat.qualifiers[tag_name].append(tag_value)
    else:
        feat.qualifiers[tag_name] = [tag_value]


def putative_lipobox(record, feature):
    tmpseq = str(feature.extract(record.seq).translate(table=11)).replace("*", "")
    return LIPOBOX.search(tmpseq)


def renumber_genes(gbk_files, problematic=None):
    for gbk_file in gbk_files:
        for record in SeqIO.parse(gbk_file, "genbank"):
            # 1. Take all the NONEs and sort out all those with genome size <90 kb.
            # 2. Call this the SMALLNONEs.  Everything else is PROBLEMATIC
            # (because lysis cassettes are not as common in large genomes)
            if len(record.seq) > 90000:
                SeqIO.write([record], problematic, 'genbank')

            # 3. Then sort out all the genomes with the words holin and
            # endolysin anywhere in the record.
            putative_endolysins = list(contains_keyword(record.features, ENDOLYSIN_KEYWORDS))

            if len(putative_endolysins) == 0:
                SeqIO.write([record], problematic, 'genbank')
            else:
                for feature in putative_endolysins:
                    tag(feature, 'color', '255 0 0')
                    tag(feature, 'note', 'putative endolysin')

                    for nearby in feat_contains_keyword(
                            features_in_range(record.features, feature, dist=2000),
                            HOLIN_KEYWORDS):
                        tag(nearby, 'color', '0 255 0')
                        tag(nearby, 'note', 'putative holin')

                    for nearby in features_in_range(record.features, feature, dist=2000):
                        if putative_lipobox(record, nearby):
                            tag(nearby, 'color', '0 0 255')
                            tag(nearby, 'note', 'putative lipobox')
                SeqIO.write([record], sys.stdout, 'genbank')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Renumber genbank files')
    parser.add_argument('gbk_files', type=file, nargs='+', help='Genbank files')
    parser.add_argument('--problematic', type=argparse.FileType('w'),
                        help='Problematic genomes')

    args = parser.parse_args()
    for record in renumber_genes(**vars(args)):
        SeqIO.write(record, sys.stdout, 'genbank')
