#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from logger import notify

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def correct_model(genbank_file):
    for record in SeqIO.parse(genbank_file, "genbank"):
        cds = [f for f in record.features if f.type == 'CDS']
        genes = [f for f in record.features if f.type == 'gene']

        # No genes at all
        if len(genes) == 0:
            genes = [
                SeqFeature(
                    f.location,
                    type="gene",
                    strand=f.location.strand,
                    qualifiers={
                        'gene': f.qualifiers.get('gene', []),
                        'locus_tag': f.qualifiers.get('locus_tag', []),
                    }
                )
                for f in cds
            ]
            # print genes
            record.features += genes
        elif len(genes) != len(cds):
            log.warn("Different number of CDSs and genes. I don't know how to handle this case (yet)")
            import StringIO
            # Eh? This is odd.
            c = StringIO.StringIO()
            SeqIO.write([record], c, 'genbank')
            notify(
                'correct_gene_model found something interesting',
                'Genes: %s, CDSs: %s' % (len(genes), len(cds)),
                {
                    'content': c.getvalue(),
                    'filename': '%s.gbk' % record.id,
                }
            )

        record.features = sorted(record.features, key=lambda x: x.location.start)
        yield [record]


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Export a subset of features from a Genbank file', epilog="")
    parser.add_argument('genbank_file', type=file, help='Genbank file')

    args = vars(parser.parse_args())
    for seq in correct_model(**args):
        SeqIO.write(seq, sys.stdout, 'genbank')
