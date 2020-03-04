#!/usr/bin/env python
import sys
import argparse
import json
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()


def parse_blast(blast):
    res = []
    finalRes = []
    for line in blast:
        taxSplit = []
        preTaxSplit = line.strip('\n').split('\t')
        for tax in preTaxSplit[-1].split(';'):
          taxSplit.append(preTaxSplit)
          taxSplit[-1][-1] = tax
          res.append(taxSplit[-1])
    for line in res:
        shallowCopy = line
        for access in line[6].split(';'):
          shallowCopy[6] = access
          finalRes.append(shallowCopy)  
    return res


def add_dice(blast):
    res = []
    for data in blast:
        dice = 2 * int(data[2]) / (float(data[3]) + float(data[4]))
        res.append(data + [dice])
    return res

    

def make_num(blast):
    res = []
    for data in blast:
        res.append([data[0], int(data[1]), int(data[2]), int(data[3]), int(data[4]), data[5], data[6], int(data[7])])
    return res

def bundle_dice(blast):
    res = []
    ind = 0
    seen = {}
    for x in blast:
      if (x[0] + x[6]) in seen.keys():
        res[seen[(x[0] + x[6])]][1] += x[1]
        res[seen[(x[0] + x[6])]][2] += x[2]
        res[seen[(x[0] + x[6])]][8] += x[8]
        res[seen[(x[0] + x[6])]][9] += 1 # Num HSPs
      else:
        seen[(x[0] + x[6])] = ind
        res.append(x + [1])
        ind += 1
    return res
        
"""def bundle_dice(blast):
    res = []
    ind = 0
    seen = {}
    for x in blast:
      if ((x[0] + x[5])) in seen.keys():
        res[seen[(x[0] + x[5])]][1] += x[1]
        res[seen[(x[0] + x[5])]][2] += x[2]
        res[seen[(x[0] + x[5])]][9] += 1 # Num HSPs
        res[seen[(x[0] + x[5])]][8] = ((res[seen[(x[0] + x[5])]][8] * (res[seen[(x[0] + x[5])]][9] - 1)) + x[8]) / res[seen[(x[0] + x[5])]][9] 
      else:
        seen[(x[0] + x[5])] = ind
        res.append(x + [1])
        ind += 1
    return res """

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
            data[0], # 01 Query Seq-id (ID of your sequence)
            data[1], # 13 All subject Seq-id(s), separated by a ';'
            split_identifiers(data[1], data[2]),  # 25 All subject title(s), separated by a '<>'
            data[3].split(";"), # Extra: All Subject Accessions
            data[4].split(";"), # Extra: All TaxIDs
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
        for ID in data[7]:
            yield [
                data[0],
                data[1],
                data[2],
                data[3],
                data[4],
                data[5],
                data[6],
                int(ID),
                data[8]
            ]

def expand_titles(blast):
    for data in blast:
        for title in data[5]:
            yield [
                data[0],
                data[1],
                data[2],
                data[3],
                data[4],
                title,
                data[6],
                data[7],
                data[8]
            ]


def filter_phage(blast, phageTaxLookup):
    res = []
    for data in blast:
        if int(data[7]) in phageTaxLookup:
            res.append(data)
    return res


def remove_dupes(data):
    has_seen = {}
    res = []
    for row in data:
        # qseqid, sseqid
        key = (row[0], row[5])
        # If we've seen the key before, we can exit
        if key in has_seen:
            continue

        # Otherwise, continue on
        has_seen[key] = True
        # Pretty simple
        res.append(row)
    return res


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

    data = [] # Reformatting to list rather than generator

    data = parse_blast(args.blast)
    data = make_num(data)
    data = add_dice(data)
    data = bundle_dice(data)
    #data = filter_dice(data, threshold=0.0)
    #data = important_only(data, splitId)
    
    #data = expand_taxIDs(data)
    #data = deform_scores(data)
    data = filter_phage(data, phageTaxLookup)
    #data = expand_titles(data)

    if args.protein or args.canonical:
        data = remove_dupes(data)  # Probably obsolete, bundle dice should do this
        count_label = "Similar Unique Proteins"
    else:
        count_label = "Nucleotide Hits"
    #data = with_dice(data)
    data.sort(key = lambda data: -data[8])
    #counts, accessions = scoreMap(data)
    
    if args.access:
      sys.stdout.write('Top %d matches for BLASTn results of %s\t\t\t\t\t\t\n' % (args.hits, data[0][0]))
      sys.stdout.write('TaxID\tName\tAccessions\tSubject Length\tNumber of HSPs\tTotal Aligned Length\tDice Score\n')
      ind = 0
      for out in data:
        if ind >= args.hits:
            break
        ind += 1
        sys.stdout.write('%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n' % (out[7], out[5], out[6], out[4], out[9], out[2], out[8]))
    else:
      sys.stdout.write('Top %d matches for BLASTn results of %s\t\t\t\t\t\n' % (args.hits, data[0][0]))
      sys.stdout.write('TaxID\tName\tSubject Length\tNumber of HSPs\tTotal Aligned Length\tDice Score\n')
      ind = 0
      for out in data:
        if ind >= args.hits:
            break
        ind += 1
        sys.stdout.write('%s\t%s\t%s\t%s\t%s\t%.4f\n' % (out[7], out[5], out[4], out[9], out[1], out[8]))
    
