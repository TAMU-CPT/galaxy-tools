#!/usr/bin/env python
import sys
import argparse
import re
import json
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()


def parse_blast(blast):
    for line in blast:
        yield line.strip('\n').split('\t')


def with_dice(blast):
    for data in blast:
        dice = 2 * int(data[14]) / (float(data[22]) + float(data[23]))
        yield data + [dice]


def filter_dice(blast, threshold=0.5):
    for data in blast:
        if data[-1] > threshold:
            yield data


def split_identifiers_nucl(_, ident):
    if '<>' in ident:
        idents = ident.split('<>')
    else:
        idents = [ident]
    return idents


def split_identifiers_prot(_, ident):
    if '<>' in ident:
        idents = ident.split('<>')
    else:
        idents = [ident]
    return [
        x[x.index('[') + 1:x.rindex(']')]
        for x in idents
        # MULTISPECIES: recombination-associated protein RdgC [Enterobacteriaceae]<>RecName: Full=Recombination-associated protein RdgC<>putative exonuclease, RdgC [Enterobacter sp. 638]
        if '[' in x and ']' in x
    ]


def split_identifiers_phage(par, ident):
    par = par.replace('lcl|', '')
    par = par[0:par.index('_prot_')]
    return [par]


def important_only(blast, split_identifiers):
    for data in blast:
        yield [
            # 01 Query Seq-id (ID of your sequence)
            # 02 Subject Seq-id (ID of the database hit)
            # 03 Percentage of identical matches
            # 04 Alignment length
            # 05 Number of mismatches
            # 06 Number of gap openings
            # 07 Start of alignment in query
            # 08 End of alignment in query
            # 09 Start of alignment in subject (database hit)
            # 10 End of alignment in subject (database hit)
            data[10], # 11 Expectation value (E-value)
            # 12 Bit score
            # 13 All subject Seq-id(s), separated by a ';'
            # 14 Raw score
            # 15 Number of identical matches
            # 16 Number of positive-scoring matches
            # 17 Total number of gaps
            # 18 Percentage of positive-scoring matches
            # 19 Query frame
            # 20 Subject frame
            # 21 Aligned part of query sequence
            # 22 Aligned part of subject sequence
            # 23 Query sequence length
            # 24 Subject sequence length
            split_identifiers(data[1], data[24]), # 25 All subject title(s), separated by a '<>'
            data[25] # 26 dice
        ]


def deform_scores(blast):
    for data in blast:
        for org in data[1]:
            yield [
                data[0],
                org,
                data[2]
            ]


def filter_phage(blast, phageNameLookup):
    for data in blast:
        if data[1] in phageNameLookup:
            yield [
                data[0],
                data[1],
                phageNameLookup[data[1]],
                data[2]
            ]


def scoreMap(blast):
    m = {}
    c = {}
    for (evalue, name, id, dice) in blast:
        if (name, id) not in m:
            m[(name, id)] = 0
            c[(name, id)] = 0

        m[(name, id)] += 1 * dice
        c[(name, id)] += 1
    return m, c


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Top related genomes')
    parser.add_argument('blast', type=argparse.FileType("r"), help='Blast 25 Column Results')
    parser.add_argument('phagedb', type=argparse.FileType("r"))
    parser.add_argument('--protein', action='store_true')
    parser.add_argument('--canonical', action='store_true')

    args = parser.parse_args()

    phageDb = json.load(args.phagedb)
    if args.protein:
        splitId = split_identifiers_prot
        phageNameLookup = {k['source'].rstrip('.'): k['id'] for k in phageDb}
    elif args.canonical:
        splitId = split_identifiers_phage
        phageNameLookup = {k['source'].rstrip('.'): k['id'] for k in phageDb}
    else:
        splitId = split_identifiers_nucl
        phageNameLookup = {k['desc'].rstrip('.'): k['id'] for k in phageDb}

    data = parse_blast(args.blast)
    data = with_dice(data)
    data = filter_dice(data, threshold=0.0)
    data = important_only(data, splitId)
    data = deform_scores(data)
    data = filter_phage(data, phageNameLookup)

    scores, counts = scoreMap(data)
    sys.stdout.write('# ID\tName\tScore\tProtein Count\n')
    for idx, ((name, pid), score) in enumerate(sorted(scores.items(), key=lambda (x, y): -y)):
        if idx > 4:
            break

        sys.stdout.write('%s\t%s\t%05.3f\t%d\n' % (pid, name, score, counts[(name, pid)]))
