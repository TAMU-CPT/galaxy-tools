#!/usr/bin/env python
import sys
import argparse
import json
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()


def parse_blast(blast):
    for line in blast:
        yield line.strip("\n").split("\t")


def with_dice(blast):
    for data in blast:
        dice = 2 * int(data[14]) / (float(data[22]) + float(data[23]))
        yield data + [dice]


def filter_dice(blast, threshold=0.5):
    for data in blast:
        if data[-1] > threshold:
            yield data


def split_identifiers_nucl(_, ident):
    if "<>" in ident:
        idents = ident.split("<>")
    else:
        idents = [ident]
    return idents


def split_identifiers_prot(_, ident):
    if "<>" in ident:
        idents = ident.split("<>")
    else:
        idents = [ident]
    return [
        x[x.index("[") + 1 : x.rindex("]")]
        for x in idents
        # MULTISPECIES: recombination-associated protein RdgC [Enterobacteriaceae]<>RecName: Full=Recombination-associated protein RdgC<>putative exonuclease, RdgC [Enterobacter sp. 638]
        if "[" in x and "]" in x
    ]


def split_identifiers_phage(par, ident):
    par = par.replace("lcl|", "")
    par = par[0 : par.index("_prot_")]
    return [par]


def important_only(blast, split_identifiers):
    for data in blast:
        yield [
            data[0],  # 01 Query Seq-id (ID of your sequence)
            data[1],  # 13 All subject Seq-id(s), separated by a ';'
            split_identifiers(
                data[1], data[2]
            ),  # 25 All subject title(s), separated by a '<>'
            data[3].split(";"),  # Extra: All Subject Accessions
            data[4].split(";"),  # Extra: All TaxIDs
        ]


def deform_scores(blast):
    for data in blast:
        for org in data[2]:
            yield [data[0], data[1], org, data[3], data[4]]


def expand_fields(blast):
    for data in blast:
        for x in range(0, len(data[4])):
            yield [data[0], data[1], data[2][x], data[3], int(data[4][x])]

def expand_taxIDs(blast):
    for data in blast:
        # if(len(data[4]) > 0):
        #  print(data[0])
        for ID in data[4]:
            if ID != "N/A":
              yield [data[0], data[1], data[2], data[3], int(ID)]


def expand_titles(blast):
    for data in blast:
        for title in data[2]:
            yield [data[0], data[1], title, data[3], data[4]]


def filter_phage(blast, phageTaxLookup, phageSciNames):
    for data in blast:
        for x in range(0, len(phageTaxLookup)):
            if (data[4]) == phageTaxLookup[x]:
                yield [data[0], data[1], phageSciNames[x], data[3], data[4]]
                break


def remove_dupes(data):
    has_seen = {}
    for row in data:
        # qseqid, sseqid
        key = (row[0], row[4])
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
        if (str(subTitle), ID) not in c:
            m[(str(subTitle), ID)] = access
            c[(str(subTitle), ID)] = 0

        c[(str(subTitle), ID)] += 1
    return c, m


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Top related genomes")
    parser.add_argument(
        "blast", type=argparse.FileType("r"), help="Blast 25 Column Results"
    )
    parser.add_argument("phagedb", type=argparse.FileType("r"))
    parser.add_argument("--access", action="store_true")
    parser.add_argument("--protein", action="store_true")
    parser.add_argument("--canonical", action="store_true")
    parser.add_argument("--noFilter", action="store_true")
    #parser.add_argument("--title", action="store_true") # Add when ready to update XML after semester
    parser.add_argument("--hits", type=int, default=5)
    

    args = parser.parse_args()

    phageDb = args.phagedb
    phageTaxLookup = []
    sciName = []
    line = phageDb.readline()
    while line:
        line = line.split("\t")
        phageTaxLookup.append(int(line[0]))
        line[1] = line[1].strip()
        if (line[1] == ""):
            line[1] = "Novel Genome"
        sciName.append(line[1])
        line = phageDb.readline()

    line = args.blast.readline()
    line = line.split("\t")
    nameRec = line[0]
    args.blast.seek(0)

    if args.protein:
        splitId = split_identifiers_prot
        # phageNameLookup = {k['source'].rstrip('.'): k['id'] for k in phageDb}
    elif args.canonical:
        splitId = split_identifiers_phage
        # phageNameLookup = {k['source'].rstrip('.'): k['id'] for k in phageDb}
    else:
        splitId = split_identifiers_nucl
        # phageNameLookup = {k['desc'].rstrip('.'): k['id'] for k in phageDb}

    data = parse_blast(args.blast)
    # data = with_dice(data)
    # data = filter_dice(data, threshold=0.0)
    data = important_only(data, splitId)
    
    data = expand_taxIDs(data)
    data = remove_dupes(data)
    if not args.noFilter:
        data = filter_phage(data, phageTaxLookup, sciName)
    listify = []
    for x in data:
        listify.append(x)
    #listify = greatest_taxID(listify)
       
    count_label = "Similar Unique Proteins"
    
    counts, accessions = scoreMap(listify)
    
    sys.stdout.write(
            "Top %d matches for BLASTp results of %s\n"
            % (args.hits, nameRec)
        )
    header = "# TaxID\t"
    #if args.title:
    header += "Name\t"
    if args.access:
        header += "Accessions\t"
    header += "Similar Unique Proteins\n"
    sys.stdout.write(header)

    for idx, ((name, ID), num) in enumerate(
            sorted(counts.items(), key=lambda item: -item[1])
        ):
        if idx > args.hits - 1:
            break
        line = str(ID) + "\t"
        #if args.title:
        line += str(name) + "\t"
        if args.access:
          line += str(accessions[(name, ID)][0]) + "\t"
        line += str(num) + "\n" 
        sys.stdout.write(line)
