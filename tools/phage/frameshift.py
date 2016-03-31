#!/usr/bin/env python
import argparse
import copy
import re
import sys
from collections import Counter
from Levenshtein import distance
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from gff3 import feature_lambda, feature_test_type


import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

def slipperyScore(sequence):
    a = Counter(str(sequence[0:3]).upper())
    b = Counter(str(sequence[3:6]).upper())
    # There are a couple of cases. 1) all three are the same, len(keys) == 1,
    # one char is diff, so len(keys) == 2, and then all are diff. Long story
    # short, can use keyset lengths here.
    score = 4 + (1 - len(a.keys())) + (1 - len(b.keys()))
    return score


def FrameShiftFinder(gff3, fasta, max_overlap=60, table=11, slippage_max=-3):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    for rec in GFF.parse(gff3, base_dict=seq_dict):
        putative_frameshift_genes = []

        for gene in feature_lambda(rec.features, feature_test_type,
                {'type': 'gene'}):
            feats = list(feature_lambda(gene.sub_features, feature_test_type, {'type': 'CDS'}))
            if len(feats) == 0:
                continue

            # only the first CDS. Bad idea...
            feat = feats[0]

            # For each feature, we check the last N bases
            overlap = min(len(feat), max_overlap)
            # Make sure overlap is evenly divisible by three
            if overlap % 3 != 0:
                overlap -= (overlap % 3)

            feat_seq = feat.extract(rec)
            # Loop across the last overlap / 3 codons
            for idx in range(len(feat) - overlap, len(feat), 3):
                codons = feat_seq[idx - 6:idx]
                # Look from -10 to +5, and check under wobble rules, one
                # or two codons are equal.
                tn_codons = codons.seq.translate(table=table)
                for wobble in range(slippage_max, 3):
                    # Don't allow wobbling within frame
                    if wobble % 3 == 0:
                        continue

                    cmp_codons = feat_seq[idx - 6 + wobble:idx+wobble]
                    cmp_codons = cmp_codons.seq.translate(table=table)

                    # Here we've found a feature which has a possible
                    # frame shift target. Now we need to check if there's
                    # a run of stop-codon free.
                    if feat.location.strand > 0:
                        possible_frameshift_start = feat.location.start + \
                            idx - 6 + wobble
                        frameshift_to_end = rec.seq[possible_frameshift_start:]
                    else:
                        possible_frameshift_start = feat.location.end - \
                            idx + 6 - wobble
                        frameshift_to_end = rec.seq[0:possible_frameshift_start].reverse_complement()

                    translated_seq = frameshift_to_end.translate(table=11, to_stop=True)

                    # Skip sequences which are too short. Divide by three
                    # as trnaslated seq is protein, and feat is DNA.
                    # Divide by 2 again so we collect sequences which are
                    # at least half as long as the parent.
                    if len(translated_seq) < len(feat) / (3 * 2):
                        continue

                    # If both AAs are different, continue. This would be very strange.
                    if distance(str(tn_codons), str(cmp_codons)) > 1:
                        continue

                    # If the edit distance between the amino acid sequences is > 4, remove it as it is unlikely.
                    if distance(str(codons.seq), str(frameshift_to_end[0:6])) > 4:
                        continue


                    # Do not allow cases where the leading NT changes.
                    # if str(codons.seq)[0] != str(cmp_codons)[0]:
                        # continue

                    # Scoring to help users filter out less likely results.
                    # Start at 4 because they've already been identified at
                    # some level as being a possible frameshift candidate.
                    score = 4 + slipperyScore(frameshift_to_end[0:7])


                    sys.stderr.write('>' + '\t'.join(map(str, (
                        score,
                        distance(str(tn_codons), str(cmp_codons)),
                        feat.id,
                        feat.location.strand,
                        possible_frameshift_start, idx,
                        wobble, codons.seq, tn_codons, cmp_codons,
                        frameshift_to_end[0:6]
                    ))) + '\n')

                    # Ok, we need to add a new mRNA structure for this
                    if feat.location.strand > 0:
                        cdsB_start = possible_frameshift_start
                        cdsB_end = possible_frameshift_start + (len(translated_seq) * 3)
                    else:
                        cdsB_end = possible_frameshift_start
                        cdsB_start = possible_frameshift_start - (len(translated_seq) * 3)

                    frameshiftedFeatureLocation = FeatureLocation(
                        cdsB_start,
                        cdsB_end,
                        strand=feat.location.strand
                    )

                    mRNA = SeqFeature(
                        frameshiftedFeatureLocation,
                        type='mRNA',
                        qualifiers={
                            'source': 'CPT_FSFinder',
                            'score': score * (1000 / 8),
                        }
                    )
                    exon = SeqFeature(
                        frameshiftedFeatureLocation,
                        type='exon',
                        qualifiers={
                            'source': 'CPT_FSFinder',
                        }
                    )
                    cds_b = SeqFeature(
                        frameshiftedFeatureLocation,
                        type='CDS',
                        qualifiers={
                            'source': 'CPT_FSFinder',
                            'phase': wobble % 3, # WHO KNOWS
                        }
                    )

                    mRNA.sub_features = [cds_b, exon]

                    gene2 = copy.deepcopy(gene)
                    # Strip the ID
                    del gene2.qualifiers['ID']
                    gene2.id = None
                    # Add more quals
                    gene2.qualifiers.update({
                        'source': 'CPT_FSFinder',
                        'Name': 'Frameshifted ' + gene2.qualifiers.get('Name', [''])[0],
                        'Note': ['Predicted frameshift region', 'Frameshift is ' + str(wobble)]
                    })
                    gene2.sub_features = [mRNA]
                    putative_frameshift_genes.append(gene2)

        rec.features = putative_frameshift_genes
        rec.annotations = {}
        yield [rec]



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gff3', type=file, help='GFF3 gene calls')
    parser.add_argument('fasta', type=file, help='FASTA genome calls')
    parser.add_argument('--max_overlap', type=int, default=60)
    parser.add_argument('--table', type=int, default=11)

    args = parser.parse_args()
    for rec in FrameShiftFinder(**vars(args)):
        GFF.write(rec, sys.stdout)
