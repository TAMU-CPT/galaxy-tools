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
            data[0], # 01 Query Seq-id (ID of your sequence)
            data[1], # 13 All subject Seq-id(s), separated by a ';'
            split_identifiers(data[1], data[2]),  # 25 All subject title(s), separated by a '<>'
            data[3].split(";"), # Extra: All Subject Accessions
            data[4].split(";") # Extra: All TaxIDs
        ]


def deform_scores(blast):
    for data in blast:
        for org in data[2]:
            yield [
                data[0],
                data[1],
                org,
                data[3],
                data[4]
            ]

def expand_taxIDs(blast):
    for data in blast:
        #if(len(data[4]) > 0):
        #  print(data[0])
        for ID in data[4]:
            yield [
                data[0],
                data[1],
                data[2],
                data[3],
                int(ID)
            ]

def expand_titles(blast):
    for data in blast:
        for title in data[2]:
            yield [
                data[0],
                data[1],
                title,
                data[3],
                data[4]
            ]


def filter_phage(blast, phageTaxLookup):
    for data in blast:
        if (data[4]) in phageTaxLookup:
            yield [
                data[0],
                data[1],
                data[2],
                data[3],
                data[4]
            ]


def remove_dupes(data):
    has_seen = {}
    for row in data:
        # qseqid, sseqid
        key = (row[0], row[2], row[4])
        # If we've seen the key before, we can exit
        if key in has_seen:
            continue

        # Otherwise, continue on
        has_seen[key] = True
        # Pretty simple
        yield row


def scoreMap(blast):
    c = {}
    m = {}
    for (qseq, subID, subTitle, access, ID) in blast:
        if (subTitle, ID) not in c:
            m[(subTitle, ID)] = access
            c[(subTitle, ID)] = 0

        c[(subTitle, ID)] += 1
    return c, m


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Top related genomes')
    parser.add_argument('blast', type=argparse.FileType("r"), help='Blast 25 Column Results')
    parser.add_argument('phagedb', type=argparse.FileType("r"))
    parser.add_argument('--access', action='store_true')
    parser.add_argument('--protein', action='store_true')
    parser.add_argument('--canonical', action='store_true')
    parser.add_argument('--hits', type = int, default = 5)

    args = parser.parse_args()

    phageDb = args.phagedb
    phageTaxLookup =[]
    line = phageDb.readline()
    while line:
      phageTaxLookup.append(int(line))
      line = phageDb.readline()
    
    if args.protein:
        splitId = split_identifiers_prot
        #phageNameLookup = {k['source'].rstrip('.'): k['id'] for k in phageDb}
    elif args.canonical:
        splitId = split_identifiers_phage
        #phageNameLookup = {k['source'].rstrip('.'): k['id'] for k in phageDb}
    else:
        splitId = split_identifiers_nucl
        #phageNameLookup = {k['desc'].rstrip('.'): k['id'] for k in phageDb}

    data = parse_blast(args.blast)
    #data = with_dice(data)
    #data = filter_dice(data, threshold=0.0)
    data = important_only(data, splitId)
    
    data = expand_taxIDs(data)
    #data = deform_scores(data)
    data = filter_phage(data, phageTaxLookup)
    data = expand_titles(data)

    if args.protein or args.canonical:
        data = remove_dupes(data)
        count_label = "Similar Unique Proteins"
    else:
        count_label = "Nucleotide Hits"

    counts, accessions = scoreMap(data)
    if args.access:
      sys.stdout.write('# TaxID\tName\tAccessions\t%s\n' % count_label)
      for idx, ((name, ID), num) in enumerate(sorted(counts.items(), key=lambda item: -item[1])):
        if idx > args.hits - 1:
            break

        sys.stdout.write('%s\t%s\t%s\t%d\n' % (ID, name, str(accessions[(name, ID)][0]), num))
    else:
      sys.stdout.write('# TaxID\tName\t%s\n' % count_label)
      for idx, ((name, ID), num) in enumerate(sorted(counts.items(), key=lambda item: -item[1])):
        if idx > args.hits - 1:
            break

        sys.stdout.write('%s\t%s\t%d\n' % (ID, name, num))
