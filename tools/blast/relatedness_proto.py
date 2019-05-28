#!/usr/bin/env python
import sys
import argparse
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
            data[0],  # 01 Query Seq-id (ID of your sequence)
            # 02 Subject Seq-id (ID of the database hit)
            # 03 Percentage of identical matches
            # 04 Alignment length
            # 05 Number of mismatches
            # 06 Number of gap openings
            # 07 Start of alignment in query
            # 08 End of alignment in query
            # 09 Start of alignment in subject (database hit)
            # 10 End of alignment in subject (database hit)
            float(data[10]),  # 11 Expectation value (E-value)
            float(data[11]),  # 12 Bit score
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
            split_identifiers(data[1], data[24]),  # 25 All subject title(s), separated by a '<>'
            data[25]  # 26 dice
        ]


def deform_scores(blast):
    for data in blast:
        for org in data[3]:
            yield [
                data[0],
                data[1],
                data[2],
                org,
                data[4]
            ]


def filter_phage(blast, phageNameLookup):
    for data in blast:
        if data[3] in phageNameLookup:
            yield [
                data[0],
                data[1],
                data[2],
                data[3],
                phageNameLookup[data[3]],
                data[4]
            ]


def remove_dupes(data):
    has_seen = {}
    for row in data:
        # qseqid, sseqid
        key = (row[0], row[3])
        # If we've seen the key before, we can exit
        if key in has_seen:
            continue

        # Otherwise, continue on
        has_seen[key] = True
        # Pretty simple
        yield row


def scoreMap(blast):
    m = {}
    c = {}
    lowE = {}
    highE = {}
    lowB = {}
    highB = {}
    for (qseq, evalue, bvalue, name, id, dice) in blast:
        if (name, id) not in m:
            m[(name, id)] = 0
            c[(name, id)] = 0
            lowE[(name, id)] = evalue
            highE[(name, id)] = evalue
            lowB[(name, id)] = bvalue
            highB[(name, id)] = bvalue

        if evalue < lowE[(name, id)]:
            lowE[(name, id)] = evalue
        elif evalue > highE[(name, id)]:
            highE[(name, id)] = evalue
        
        if bvalue < lowB[(name, id)]:
            lowB[(name, id)] = bvalue
        elif bvalue > highB[(name, id)]:
            highB[(name, id)] = bvalue

        m[(name, id)] += 1 * dice
        c[(name, id)] += 1
    return m, c, lowE, highE, lowB, highB


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Top related genomes')
    parser.add_argument('blast', type=argparse.FileType("r"), help='Blast 25 Column Results')
    parser.add_argument('phagedb', type=argparse.FileType("r"))
    parser.add_argument('--protein', action='store_true')
    parser.add_argument('--canonical', action='store_true')
    parser.add_argument('--hits', type = int, default = 5)

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
    if args.protein or args.canonical:
        data = remove_dupes(data)
        count_label = "Similar Unique Proteins"
    else:
        count_label = "Nucleotide Hits"

    scores, counts, lowEs, highEs, lowBs, highBs = scoreMap(data)
    sys.stdout.write('# ID\tName\tScore\t%s\tLowest E-Value\tHighest E-Value\tLowest Bit Value\tHighest Bit Value\n' % count_label)
    for idx, ((name, pid), score) in enumerate(sorted(scores.items(), key=lambda (x, y): -y)):
        if idx > args.hits - 1:
            break

        sys.stdout.write('%s\t%s\t%05.3f\t%d\t%5.25f\t%5.25f\t%5.5f\t%5.5f\n' % (pid, name, score, counts[(name, pid)], lowEs[(name, pid)], highEs[(name, pid)], lowBs[(name, pid)], highBs[(name, pid)]))
