#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import numpy
import logging
from Bio import Entrez
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


__doc__ = """
Find top related genomes
========================

Phage genomes often share large percentages of their genome. This tool allows
easy determination of related genomes from blast results (must be 12+ column tabular)
"""


def reduce_to_score(evalue_list):
    # Smaller than 1e-250 just goes to zero afaict
    raw = [numpy.abs(numpy.log(x+1e-250)) for x in evalue_list]
    # Compensate for length somewhat
    score = numpy.sum(raw)# * len(raw)
    return (score, {
        'num': len(raw),
        'mean': numpy.mean(raw),
        'median': numpy.median(raw),
        'std': numpy.std(raw),
    })


def top_related(blast=None, method='n', restrict_to_phage=True, email=None, **kwd):
    # Table of hits
    hits = {}
    # Keep a table of example GI numbers, so we can hopefully do reverse genome lookups from the protein
    giex = {}
    for line in blast.readlines():
        split_line = line.split('\t')

        # Important data
        evalue = float(split_line[10])
        gi_hits = split_line[12].split(';')
        org_hits = []
        for x in [hit for hit in split_line[24].strip('MULTISPECIES: ').split('<>')]:
            try:
                # No regex here because of http://www.ncbi.nlm.nih.gov/protein/355386069
                #
                # recombination protein U [ [[Clostridium] clostridioforme 2_1_49FAA]]
                #
                # Which appears precisely as above in the blast output.
                org_hits.append(hit[hit.index('[')+1:hit.rindex(']')].strip())
            except:
                org_hits.append(hit.strip())

        # Thanks to Peter's parser, the gi list and org list are the same
        # length (the first gi column is also the first gi in the "master" gi
        # column)
        for (org, gi) in zip(org_hits, gi_hits):
            if (restrict_to_phage and ('phage' in org or 'virus' in org or 'Phage' in org or 'Bacteriophage' in org)) or not restrict_to_phage:
                if org in hits:
                    hits[org].append(evalue)
                else:
                    # For new hits store evalues and use as an example GI for that
                    # organism
                    hits[org] = [evalue]
                    giex[org] = gi

    # Reduce to an easily sortable score
    extra_data = {}
    for item in hits:
        (score, extra) = reduce_to_score(hits[item])
        hits[item] = score
        extra_data[item] = extra


    # Top results
    top_accessions = {}
    # Allow different result process methods
    if method == 'sd':
        std = numpy.std(hits.values())
        mean = numpy.mean(hits.values())
        for key, value in hits.iteritems():
            if value > (mean + 2 * std):
                top_accessions[key] = value
    else:
        for key, value in reversed(sorted(hits.iteritems(), key=lambda (k, v):
                                          (v, k))):

            if len(top_accessions) > 10:
                break
            else:
                top_accessions[key] = value

    top_names = {
        'Sheet1': {
                'header': ['Organism', 'Score', 'Number of hits', 'Mean', 'Median', 'Std Dev', 'Genome from same Taxonomy'],
                'data': [],
        }
    }

    acc_list = []
    for hit in top_accessions:
        acc = get_refseq_for_name(name=hit, email=email)
        acc_list.append(acc)
        top_names['Sheet1']['data'].append([
            hit,
            hits[hit],
            extra_data[hit]['num'],
            extra_data[hit]['mean'],
            extra_data[hit]['median'],
            extra_data[hit]['std'],
            acc,
        ])
    return (acc_list, top_names)


def get_refseq_for_name(name=None, email=None):
    Entrez.email = email
    # Find entries matching the query
    searchResultHandle = Entrez.esearch(db='nuccore', dbfrom='genome', term='%s[Organism] and ("complete genome" or "complete sequence")' % name)
    searchResult = Entrez.read(searchResultHandle)
    searchResultHandle.close()
    results = searchResult['IdList']

    # If no results, drop complete genome/complete sequence
    if len(results) == 0:
        searchResultHandle = Entrez.esearch(db='nuccore', dbfrom='genome', term='%s[Organism]' % name)
        searchResult = Entrez.read(searchResultHandle)
        searchResultHandle.close()
        results = searchResult['IdList']

    if len(searchResult['IdList']) > 1:
        log.warn("Found %s results for %s" % (len(searchResult['IdList']), name))

    if len(searchResult['IdList']) > 0:
        return searchResult['IdList'][0]
    else:
        log.warn("No results found")
        return []



if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['blast', 'Blast results', {'required': True, 'validate':
                                        'File/Input'}],
            ['restrict_to_phage', 'Restrict results to Phage/Viruses', {'validate': 'Flag'}],
            ['method', 'Method to select "top" results',
             {'required': True, 'validate': 'Option', 'default': 'n',
              'options': {'n': 'Top 5 results', 'sd':
                          'Statistical significance'}}],
            ['email', 'Email', {'required': True, 'default': 'cpt@tamu.edu', 'validate': 'String'}],
        ],
        outputs=[
            [
                'accession_list',
                'Accession numbers of top matched genomes',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'top_accessions',
                    'data_format': 'text/plain',
                    'default_format': 'TXT',
                }
            ],
            [
                'name_list',
                'Human readable names of matched genomes',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'top_names',
                    'data_format': 'text/tabular',
                    'default_format': 'TSV_U',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.genbank.TopRelated',
            'appname': 'Top Related Genomes',
            'appvers': '0.7',
            'appdesc': 'finds top related sequences from blast tsv results',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    (txt, names) = top_related(**options)

    from galaxygetopt.outputfiles import OutputFiles
    of2 = OutputFiles(name='accession_list', GGO=opts)
    of2.CRR(data='\n'.join(txt))
    of3 = OutputFiles(name='name_list', GGO=opts)
    of3.CRR(data=names)
