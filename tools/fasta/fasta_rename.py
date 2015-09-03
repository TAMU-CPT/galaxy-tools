#!/usr/bin/env python
import logging
logging.basicConfig(level=logging.INFO)
import argparse


def rename_fasta_sequences(fasta_file, new_name):
    from Bio import SeqIO
    import StringIO
    output = StringIO.StringIO()

    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) > 1:
        raise Exception("Too many sequences")
    elif len(records) == 0:
        raise Exception("Too few sequences")

    orig = records[0].id
    records[0].id = new_name
    records[0].description = " [Orig=%s]" % orig
    fasta_file.close()
    SeqIO.write(records, output, "fasta")
    return output.getvalue()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rename fasta sequences')
    parser.add_argument('fasta_file', metavar='N', type=file, nargs='?',
                        help='fasta file')
    parser.add_argument('new_name', nargs='?', help='New name for the fasta sequence')
    args = parser.parse_args()

    print rename_fasta_sequences(**vars(args))
