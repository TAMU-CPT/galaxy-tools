#!/usr/bin/env python
import numpy
import argparse
import re
from kyotocabinet import DB
import logging
logging.basicConfig(level=logging.DEBUG)
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
    # Connect to kyoto db
    db = DB()
    if not db.open("/opt/gene2accession/gene2accession.kch", DB.OWRITER | DB.OCREATE):
        raise Exception("Could not load gene2accession.kch: " + str(db.error()))

    hits = {}
    gi_num = re.compile('gi\|([0-9]+)')
    for line in blast:
        split_line = line.split('\t')

        # Important data
        evalue = float(split_line[10])

        gi_nums = gi_num.findall(split_line[12])
        genome_ids = [db.get(x) for x in gi_nums if db.get(x) is not None]

        # Thanks to Peter's parser, the gi list and org list are the same
        # length (the first gi column is also the first gi in the "master" gi
        # column)
        for org in genome_ids:
            if org in hits:
                hits[org].append(evalue)
            else:
                hits[org] = [evalue]
    db.close()
    return hits


def top_related(blast, report=None):
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

    print "# " + '\t'.join(['Accession', 'Score', 'Number of hits', 'Mean',
                            'Median', 'Std Dev'])

    for hit in top_accessions:
        print '\t'.join(map(str, [
            hit,
            hits[hit],
            extra_data[hit]['num'],
            extra_data[hit]['mean'],
            extra_data[hit]['median'],
            extra_data[hit]['std'],
        ]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Top related genomes')
    parser.add_argument('blast', type=file, help='Blast 25 Column Results')

    args = parser.parse_args()
    top_related(**vars(args))
