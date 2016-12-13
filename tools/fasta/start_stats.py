#!/usr/bin/env python
import argparse
from Bio import SeqIO


def extract_starts(fasta):
    codon_usage = {}
    for record in SeqIO.parse(fasta, 'fasta'):
        seq = record.seq[0:3]
        sseq = str(seq)
        try:
            codon_usage[sseq] = (codon_usage[sseq][0] + 1, seq)
        except KeyError:
            codon_usage[sseq] = (1, seq)

    for (sseq, (count, seq)) in sorted(codon_usage.items()):
        yield (
            sseq,
            codon_usage[sseq][0],
            seq.translate(table=11)
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarise start codon usage', epilog="")
    parser.add_argument('fasta', type=argparse.FileType("r"), help='Fasta Genome')
    args = parser.parse_args()

    print('# DNA\tCodon\tCount')
    for (key, value, codon) in extract_starts(**vars(args)):
        print('{}\t{}\t{}'.format(key, codon, value))
