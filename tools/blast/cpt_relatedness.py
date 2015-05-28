#!/usr/bin/env python
import numpy
import argparse
from Bio import Entrez

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def reduce_to_score(evalue_list):
    # Smaller than 1e-250 just goes to zero afaict
    raw = [numpy.abs(numpy.log(x + 1e-250)) for x in evalue_list]
    # Compensate for length somewhat
    score = numpy.sum(raw)  # * len(raw)
    return (score, {
        'num': len(raw),
        'mean': numpy.mean(raw),
        'median': numpy.median(raw),
        'std': numpy.std(raw),
    })


def __load_blast_data(blast):
    hits = {}
    for line in blast:
        split_line = line.split('\t')

        # Important data
        evalue = float(split_line[10])
        org_hits = []
        for x in [hit for hit in split_line[24].strip('MULTISPECIES: ').split('<>')]:
            if '[' in hit and ']' in hit:
                # recombination protein U [ [[Clostridium] clostridioforme 2_1_49FAA]]
                #
                # Use rindex of [ and index of ] instead of vice versa
                # so the above string will be picked up as "Clostridium"
                org_hits.append(hit[hit.rindex('[') + 1:hit.index(']')].strip())
            # Excluding things without proper "protein name [organism name]"
            # strings improves search result quality **drastically**
            #
            # org_hits.append(hit.strip())

        # Thanks to Peter's parser, the gi list and org list are the same
        # length (the first gi column is also the first gi in the "master" gi
        # column)
        for org in org_hits:
            if org in hits:
                hits[org].append(evalue)
            else:
                hits[org] = [evalue]

    return hits


def top_related(blast, email, report=None):
    # hits = Table of hits
    hits = __load_blast_data(blast)

    # Reduce to an easily sortable score
    extra_data = {}
    for item in hits:
        (score, extra) = reduce_to_score(hits[item])
        hits[item] = score
        extra_data[item] = extra

    # Top results
    top_accessions = {}
    for key, value in list(reversed(sorted(hits.iteritems(), key=lambda (k, v):
                                           (v, k))))[0:10]:

        top_accessions[key] = value

    report.write('\t'.join(['Organism', 'Score', 'Number of hits', 'Mean',
                            'Median', 'Std Dev', 'Genome from same Taxonomy']) + "\n")

    acc_list = []
    for hit in top_accessions:
        acc = get_refseq_for_name(name=hit, email=email)
        if isinstance(acc, str):
            acc_list.append(acc)

        report.write('\t'.join(map(str, [
            hit,
            hits[hit],
            extra_data[hit]['num'],
            extra_data[hit]['mean'],
            extra_data[hit]['median'],
            extra_data[hit]['std'],
            acc,
        ])) + '\n')

    return acc_list


def get_refseq_for_name(name=None, email=None):
    # Find entries matching the query
    completes = ('complete genome', 'complete sequence', 'whole genome', 'whole sequence')
    query = '"%s"[ORGN] AND (%s) AND gbdiv phg [PROP]' % (name, ' OR '.join(['"%s"[TITLE]' % x for x in completes]))
    log.info(query)

    searchResultHandle = Entrez.esearch(db='nuccore', dbfrom='genome', term=query)
    searchResult = Entrez.read(searchResultHandle)
    searchResultHandle.close()
    results = searchResult['IdList']

    # If no results, drop complete genome/complete sequence
    if len(results) == 0:
        query = '"%s"[ORGN] AND (gbdiv phg [PROP] OR gbdiv bct [PROP])' % name
        searchResultHandle = Entrez.esearch(db='nuccore', dbfrom='genome', term=query)
        searchResult = Entrez.read(searchResultHandle)
        searchResultHandle.close()
        log.info(query)
        log.warn("Failed over to secondary query")
        results = searchResult['IdList']

    if len(searchResult['IdList']) > 1:
        log.warn("Found %s results for %s" % (len(searchResult['IdList']), name))

    if len(searchResult['IdList']) > 0:
        return searchResult['IdList'][0]
    else:
        log.warn("No results found")
        return []


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Top related genomes')
    parser.add_argument('blast', type=file, help='Blast results')
    parser.add_argument('--email', help='Email for NBCI records')

    parser.add_argument('--report', type=argparse.FileType('w'),
                        help='Location to store report', default='top_related.tsv')

    args = parser.parse_args()
    Entrez.email = args.email
    print '\n'.join(top_related(**vars(args)))
