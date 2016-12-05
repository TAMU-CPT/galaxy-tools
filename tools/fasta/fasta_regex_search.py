#!/usr/bin/env python
import sys
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def search(fasta, regex, stats_file=None):
    if stats_file is not None:
        stats_file.write('\t'.join(('# Record ID', 'Hit #', 'Sequence Length', 'Hit start', 'Hit end')) + '\n')

    m = re.compile(regex)
    for record in SeqIO.parse(fasta, "fasta"):
        if m.search(str(record.seq).upper()):
            oseq = str(record.seq).upper()
            newseq = ""
            last_idx = 0
            for idx, hit in enumerate(m.finditer(oseq)):
                newseq += oseq[last_idx:hit.start()].lower()
                newseq += oseq[hit.start():hit.end()].upper()
                last_idx = hit.end()

                if stats_file is not None:
                    stats_file.write('\t'.join(map(str, (
                        record.id,
                        idx,
                        len(oseq),
                        hit.start(),
                        hit.end()
                    ))) + '\n')

            newseq += oseq[last_idx:].lower()
            # Update sequence
            record.seq = Seq(newseq)
            yield [record]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter sequences based on regular expressions')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='Fasta Sequence')
    parser.add_argument('regex', help='Regular Expression (use upper case sequences)')
    parser.add_argument('stats_file', type=argparse.FileType('w'), help='Output statistics')

    args = parser.parse_args()

    for record in search(**vars(args)):
        SeqIO.write(record, sys.stdout, "fasta")
