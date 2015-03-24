#!/usr/bin/env python
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name=__name__)
import argparse

def calculate_gc_skew(seq):
    counts = {
        'G': 0,
        'C': 0,
        'A': 0,
        'T': 0,
    }
    for char in seq.upper():
        try:
            counts[char] += 1
        except KeyError:
            log.warning("Unknown character %s" % char)

    gpc = counts['G'] + counts['C']
    gmc = counts['G'] - counts['C']
    if gpc > 0:
        return float(gmc) / float(gpc)
    return 0

def gc_wig(fasta_file, window=100, step_size=50):
    from Bio import SeqIO
    for record in SeqIO.parse(fasta_file, "fasta"):
        seqstr = str(record.seq)

        modseqstr = seqstr[-window:] + seqstr + seqstr[0:window]
        halfwindow = window / 2
        for i in range(1, len(seqstr), step_size):
            start = i + window - halfwindow
            end = i + window + halfwindow
            skew = calculate_gc_skew(modseqstr[start:end])
            print '\t'.join(map(str, [
                record.id,
                i-step_size/2 if i != 1 else i,
                i+step_size/2 if (i+step_size/2) < len(seqstr) else i,
                skew
            ]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate GC Skew', epilog="")
    parser.add_argument('fasta_file', type=file, help='fasta file')
    parser.add_argument('--window', type=int, help='Window size for skew calculation', default=100)
    parser.add_argument('--step_size', type=int, help='Step size for skew calculation', default=100)
    args = parser.parse_args()
    gc_wig(**vars(args))
