#!/usr/bin/env python
import logging
logging.basicConfig(level=logging.INFO)
import argparse
from Bio import SeqIO, Seq


def codon_stats(fasta_file, **kwargs):
    records = list(SeqIO.parse(fasta_file, "fasta"))

    codons = {}

    for b1 in ('A', 'C', 'T', 'G'):
        for b2 in ('A', 'C', 'T', 'G'):
            for b3 in ('A', 'C', 'T', 'G'):
                codons[b1+b2+b3] = str(Seq.Seq(b1+b2+b3).translate(table=1))

    tn_table_keys = sorted(codons.keys())

    header = ['#ID', 'Length'] + ['%s (%s)' % (x, codons[x]) for x in tn_table_keys]
    yield header

    for record in records:
        seq = str(record.seq)
        codon_counts = {}

        for tri_nt in [seq[i:i+3] for i in range(0, len(seq), 3)]:
            try:
                codon_counts[tri_nt] += 1
            except:
                codon_counts[tri_nt] = 1

        row = [record.id, len(record.seq)]
        numbers = []
        for tri_nt in tn_table_keys:
            if tri_nt in codon_counts:
                numbers.append(codon_counts[tri_nt])
            else:
                numbers.append(0)

        numbers = [float(x)/sum(numbers) for x in numbers]
        yield row + numbers


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate AA frequencies in sequences')
    parser.add_argument('fasta_file', type=file, help='Fasta file')
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()

    for row in codon_stats(**vars(args)):
        print '\t'.join(map(str, row))
