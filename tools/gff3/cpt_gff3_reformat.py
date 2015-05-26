#!/usr/bin/env python
import logging
logging.basicConfig(level=logging.INFO)
import argparse
from BCBio import GFF

def reformat(data):
    import StringIO
    output = StringIO.StringIO()

    recs = list(GFF.parse(data))
    GFF.write(recs, output)

    print output.getvalue()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reformat GFF files')
    parser.add_argument('data', type=file, help='Input annotations')
    args = parser.parse_args()
    reformat(**vars(args))
