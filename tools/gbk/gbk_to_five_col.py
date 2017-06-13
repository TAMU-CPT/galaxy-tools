#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import argparse
import logging
from Bio import SeqIO
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

# Read in Genbank file and parse features
# Output features into Five Column format

"""
>Feature SeqID
Line 1
    Column 1: Start location (first nucleotide) of a feature
    Column 2: Stop location (last nucleotide) of a feature
    Column 3: Feature name (for example, 'CDS' or 'mRNA' or 'rRNA' or 'gene' or 'exon')
Line2:
    Column 4: Qualifier name (for example, 'product' or 'number' or 'gene' or 'note')
    Column 5: Qualifier value

Repeat for each feature in a seq
Repeat Line 2 for each qualifier in a feature
"""

def gbk_to_5col(genbank):
    """Converts genbank to BankIt five column format"""
    for record in SeqIO.parse(genbank, 'genbank'):
        print('>Feature %s' % record.id)
        for feature in record.features:
            for index, part in enumerate(feature.location.parts):
                if part.strand > 0:
                    start = int(part.start)
                    end = int(part.end)
                else:
                    start = int(part.end)
                    end = int(part.start)
                if index == 0:
                    name = feature.type
                    print("%d\t%d\t%s" % (start, end, name))
                else:
                    print("%d\t%d" % (start, end))
            for (qualifier, values) in feature.qualifiers.items():
                for value in values:
                    print("\t\t\t%s\t%s" % (qualifier, value))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a Genbank file into five column format')
    parser.add_argument('genbank', type=argparse.FileType("r"), help='Genbank file')

    args = vars(parser.parse_args())
    gbk_to_5col(**args)
