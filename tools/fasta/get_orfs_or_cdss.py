#!/usr/bin/env python
import sys
import re
import argparse
import logging
logging.basicConfig()
log = logging.getLogger()
from Bio.Seq import Seq, reverse_complement, translate
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Data import CodonTable

class OrfFinder(object):

    def __init__(self, fasta_file, seq_format, table, ftype, ends, mode,
                 min_len, out_nuc, out_prot, out_bed, out_gff3):

        self.table = table
        self.table_obj = CodonTable.ambiguous_generic_by_id[table]
        self.ends = ends
        self.ftype = ftype
        self.min_len = min_len

        if seq_format.lower() == "sff":
            seq_format = "sff-trim"
        elif seq_format.lower() == "fasta":
            seq_format = "fasta"
        elif seq_format.lower().startswith("fastq"):
            seq_format = "fastq"

        log.debug("Genetic code table %i" % table)
        log.debug("Minimum length %i aa" % self.min_len)

        self.starts = sorted(self.table_obj.start_codons)
        self.stops = sorted(self.table_obj.stop_codons)
        self.re_starts = re.compile("|".join(self.starts))
        self.re_stops = re.compile("|".join(self.stops))

        if mode == "all":
            get_peptides = self.get_all_peptides
        elif mode == "top":
            get_peptides = self.get_top_peptides
        elif mode == "one":
            get_peptides = self.get_one_peptide

        out_count = 0

        out_gff3.write('##gff-version 3\n')

        for idx, record in enumerate(SeqIO.parse(fasta_file, seq_format)):
            for i, (f_start, f_end, f_strand, n, t) in enumerate(get_peptides(str(record.seq).upper())):
                out_count += 1

                descr = "length %i aa, %i bp, from %s..%s[%s] of %s" \
                        % (len(t), len(n), f_start, f_end, f_strand, record.description)
                fid = record.id + "|%s%i" % (self.ftype, i + 1)

                r = SeqRecord(Seq(n), id=fid, name="", description=descr)
                t = SeqRecord(Seq(t), id=fid, name="", description=descr)

                SeqIO.write(r, out_nuc, "fasta")
                SeqIO.write(t, out_prot, "fasta")

                nice_strand = '+' if f_strand == +1 else '-'

                out_bed.write('\t'.join(map(str, [
                    record.id, f_start, f_end, fid, 0, nice_strand])) + '\n')

                out_gff3.write('\t'.join(map(str, [
                    record.id, 'getOrfsOrCds', 'CDS', f_start + 1, f_end, '.',
                    nice_strand, 0, 'ID=%s%s' % (self.ftype, i + 1)])) + '\n')

        log.info("Found %i %ss in %i sequences", out_count, self.ftype, idx + 1)

    def start_chop_and_trans(self, s, strict=True):
        """Returns offset, trimmed nuc, protein."""
        if strict:
            assert s[-3:] in self.stops, s
        assert len(s) % 3 == 0
        for match in self.re_starts.finditer(s):
            # Must check the start is in frame
            start = match.start()
            if start % 3 == 0:
                n = s[start:]
                assert len(n) % 3 == 0, "%s is len %i" % (n, len(n))
                if strict:
                    t = translate(n, self.table, cds=True)
                else:
                    # Use when missing stop codon,
                    t = "M" + translate(n[3:], self.table, to_stop=True)
                return start, n, t
        return None, None, None

    def break_up_frame(self, s):
        """Returns offset, nuc, protein."""
        start = 0
        for match in self.re_stops.finditer(s):
            index = match.start() + 3
            if index % 3 != 0:
                continue
            n = s[start:index]
            if self.ftype == "CDS":
                offset, n, t = self.start_chop_and_trans(n)
            else:
                offset = 0
                t = translate(n, self.table, to_stop=True)
            if n and len(t) >= self.min_len:
                yield start + offset, n, t
            start = index
        if self.ends == "open":
            # No stop codon, Biopython's strict CDS translate will fail
            n = s[start:]
            # Ensure we have whole codons
            # TODO - Try appending N instead?
            # TODO - Do the next four lines more elegantly
            if len(n) % 3:
                n = n[:-1]
            if len(n) % 3:
                n = n[:-1]
            if self.ftype == "CDS":
                offset, n, t = self.start_chop_and_trans(n, strict=False)
            else:
                offset = 0
                t = translate(n, self.table, to_stop=True)
            if n and len(t) >= self.min_len:
                yield start + offset, n, t

    def get_all_peptides(self, nuc_seq):
        """Returns start, end, strand, nucleotides, protein.

        Co-ordinates are Python style zero-based.
        """
        # TODO - Refactor to use a generator function (in start order)
        # rather than making a list and sorting?
        answer = []
        full_len = len(nuc_seq)
        for frame in range(0, 3):
            for offset, n, t in self.break_up_frame(nuc_seq[frame:]):
                start = frame + offset  # zero based
                answer.append((start, start + len(n), +1, n, t))
        rc = reverse_complement(nuc_seq)
        for frame in range(0, 3):
            for offset, n, t in self.break_up_frame(rc[frame:]):
                start = full_len - frame - offset  # zero based
                answer.append((start - len(n), start, -1, n, t))
        answer.sort()
        return answer

    def get_top_peptides(self, nuc_seq):
        """Returns all peptides of max length."""
        values = list(get_all_peptides(nuc_seq))
        if not values:
            raise StopIteration
        max_len = max(len(x[-1]) for x in values)
        for x in values:
            if len(x[-1]) == max_len:
                yield x

    def get_one_peptide(self, nuc_seq):
        """Returns first (left most) peptide with max length."""
        values = list(get_top_peptides(nuc_seq))
        if not values:
            raise StopIteration
        yield values[0]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get open reading frames')
    parser.add_argument('fasta_file', type=file, help='Fasta file')

    parser.add_argument('-f', '--format', dest='seq_format',
                    default='fasta', help='Sequence format (e.g. fasta, fastq, sff)')
    parser.add_argument('--table', dest='table',
                    default=1, help='NCBI Translation table', type=int)
    parser.add_argument('-t', '--ftype', dest='ftype',
                    choices=('CDS', 'ORF'), default='ORF',
                    help='Find ORF or CDSs')
    parser.add_argument('-e', '--ends', dest='ends',
                    choices=('open', 'closed'), default='closed',
                    help='Open or closed. Closed ensures start/stop codons are present')
    parser.add_argument('-m', '--mode', dest='mode',
                    choices=('all', 'top', 'one'), default='all',
                    help='Output all ORFs/CDSs from sequence, all ORFs/CDSs '
                    'with max length, or first with maximum length')
    parser.add_argument('--min_len', dest='min_len',
                    default=10, help='Minimum ORF/CDS length', type=int)

    parser.add_argument('--on', dest='out_nuc',type=argparse.FileType('w'),
                    default='out.fna', help='Output nucleotide sequences')
    parser.add_argument('--op', dest='out_prot',type=argparse.FileType('w'),
                    default='out.fa', help='Output protein sequences',)
    parser.add_argument('--ob', dest='out_bed',type=argparse.FileType('w'),
                    default='out.bed', help='Output BED file')
    parser.add_argument('--og', dest='out_gff3', type=argparse.FileType('w'),
                    default='out.gff3', help='Output GFF3 file')
    parser.add_argument('-v', action='version', version='0.3.0')
    args = parser.parse_args()

    OrfFinder(**vars(args))
