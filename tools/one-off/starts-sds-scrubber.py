#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
from Bio.Data import CodonTable
from Bio.Alphabet import IUPAC
import itertools
import sys

# get all possible codons for every amino acid
table = CodonTable.unambiguous_dna_by_id[11]
aa_codes = {}
for codon in table.forward_table:
    if table.forward_table[codon] in aa_codes:
        aa_codes[table.forward_table[codon]].append(codon)
    else:
        aa_codes[table.forward_table[codon]] = [codon]

sds = (  # all possible SD sequences
    'AGGAGGT',
    'GGAGGT',
    'AGGAGG',
    'GGGGGG',
    'AGGAG',
    'GAGGT',
    'GGAGG',
    'GGGGG',
    'AGGT',
    'GGGT',
    'GAGG',
    'GGGG',
    'AGGA',
    'GGAG',
    'GGA',
    'GAG',
    'AGG',
    'GGT',
    'GGG',
)

def get_CDS_and_SD(feat):
    """
        return the CDS and SD feature in gene
    """
    for sf in feat.sub_features:
        if sf.type == 'Shine_Dalgarno_sequence':
            sd = sf
        elif sf.type == 'mRNA':
            for sfmrna in sf.sub_features:
                if sfmrna.type == 'CDS':
                    cds = sfmrna

    return cds, sd

def break_start(seq):
    """
        if possible, change the start sequence
        while keeping the same amino acid translation
    """
    seq = seq.upper()
    if seq == 'ATG':    # if methionine, no other option but ATG
        return 'ATG'
    elif seq == 'GTG':  # if valine, return GTA (still valine)
        return 'GTA'
    elif seq == 'TTG':  # if leucine, return TTA (still leucine)
        return 'TTA'

def changed_letters(a, b):
    return sum(a[i] != b[i] for i in range(len(a)))

def break_sd(sd):
    """
        if possible, change the SD sequence
        while keeping the same amino acid translation
    """
    cdns = []
    for c in sd.translate(table=table):  # find all possible combinations of codons w/ same translation
        cdns.append(aa_codes[c])

    poss_changes = []  # changes that keep same aa sequence
    for i in list(itertools.product(*cdns)):  # if a combo has no shine, return
        check = ''.join(i)
        no_sd = True
        for s in sds:
            if s in check:
                no_sd = False
        if no_sd:
            poss_changes.append(check)

    num_changes = 100
    best_replacement = sd
    for p in poss_changes:  # return the aa seq with fewest changes
        if changed_letters(p, str(sd)) < num_changes:
            num_changes = changed_letters(p, str(sd))
            best_replacement = Seq(p, IUPAC.unambiguous_dna)

    return best_replacement

def next_first_frame(start, mod_pos):
    """ return position of next nucleotide in frame 1 """
    mod = mod_pos
    while ((start - mod) % 3) + 1 != 1:
        mod += mod_pos
    return mod

def repair(fasta, gff3):
    recs = {}
    # seqids = {}
    for record in GFF.parse(gff3):
        # seqids[record.id] = ''
        recs[record.id] = record

    seqs = []
    for seq in SeqIO.parse(fasta, 'fasta'):
        if seq.id not in recs:
            continue

        current = recs[seq.id]
        for num, feat in enumerate(current.features):
            if num == 0:  # ignore first feature bc that's the full one
                continue

            cds, sd = get_CDS_and_SD(feat)
            cds_start = seq.seq[cds.location.start:cds.location.start+3]
            broken_start = break_start(cds_start)

            if cds_start != broken_start:  # try to break start sequence while keeping amino acid the same
                seq.seq = seq.seq[0:cds.location.start] + broken_start + seq.seq[cds.location.start+3:]

            else:  # if couldn't change start, must break SD
                mod_sd_start = 0
                mod_sd_end = 0
                if (sd.location.start % 3) + 1 != 1:
                    mod_sd_start = next_first_frame(sd.location.start, 1)
                if (sd.location.end % 3) + 1 != 1:
                    mod_sd_end = next_first_frame(sd.location.end, -1)

                sd_seq = seq.seq[sd.location.start-mod_sd_start:sd.location.end-mod_sd_end]
                broken_sd = break_sd(sd_seq)
                if sd_seq != broken_sd:
                    seq.seq = seq.seq[0:(sd.location.start-mod_sd_start)] + broken_sd + seq.seq[(sd.location.end-mod_sd_end):]

        seqs.append(seq)

    SeqIO.write(seqs, sys.stdout, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='change bases o get rid of internal starts/SDs')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='fasta of gene seqs')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='gff3 of starts in each seq')
    args = parser.parse_args()

    repair(args.fasta, args.gff3)
