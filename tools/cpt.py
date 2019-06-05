#!/usr/bin/env python
import regex as re
from Bio.Seq import Seq, reverse_complement, translate
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Data import CodonTable
import logging
logging.basicConfig()
log = logging.getLogger()

PHAGE_IN_MIDDLE = re.compile('^(?P<host>.*)\s*phage (?P<phage>.*)$')
BACTERIOPHAGE_IN_MIDDLE = re.compile('^(?P<host>.*)\s*bacteriophage (?P<phage>.*)$')
STARTS_WITH_PHAGE = re.compile('^(bacterio|vibrio|Bacterio|Vibrio|)?[Pp]hage (?P<phage>.*)$')
NEW_STYLE_NAMES = re.compile('(?P<phage>v[A-Z]_[A-Z][a-z]{2}_.*)')


def phage_name_parser(name):
    host = None
    phage = None
    name = name.replace(', complete genome.', '')
    name = name.replace(', complete genome', '')

    m = BACTERIOPHAGE_IN_MIDDLE.match(name)
    if m:
        host = m.group('host')
        phage = m.group('phage')
        return (host, phage)

    m = PHAGE_IN_MIDDLE.match(name)
    if m:
        host = m.group('host')
        phage = m.group('phage')
        return (host, phage)

    m = STARTS_WITH_PHAGE.match(name)
    if m:
        phage = m.group('phage')
        return (host, phage)

    m = NEW_STYLE_NAMES.match(name)
    if m:
        phage = m.group('phage')
        return (host, phage)

    return (host, phage)


class OrfFinder(object):

    def __init__(self, table, ftype, ends, min_len, strand):
        self.table = table
        self.table_obj = CodonTable.ambiguous_generic_by_id[table]
        self.ends = ends
        self.ftype = ftype
        self.min_len = min_len
        self.starts = sorted(self.table_obj.start_codons)
        self.stops = sorted(self.table_obj.stop_codons)
        self.re_starts = re.compile("|".join(self.starts))
        self.re_stops = re.compile("|".join(self.stops))
        self.strand = strand

    def locate(self, fasta_file, out_nuc, out_prot, out_bed, out_gff3):
        seq_format = "fasta"
        log.debug("Genetic code table %i" % self.table)
        log.debug("Minimum length %i aa" % self.min_len)

        out_count = 0

        out_gff3.write('##gff-version 3\n')

        for idx, record in enumerate(SeqIO.parse(fasta_file, seq_format)):
            for i, (f_start, f_end, f_strand, n, t) in enumerate(self.get_all_peptides(str(record.seq).upper())):
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
                    nice_strand, 0, 'ID=%s.%s.%s' % (self.ftype, idx, i + 1)])) + '\n')
        log.info("Found %i %ss", out_count, self.ftype)

    def start_chop_and_trans(self, s, strict=True):
        """Returns offset, trimmed nuc, protein."""
        if strict:
            assert s[-3:] in self.stops, s
        assert len(s) % 3 == 0
        for match in self.re_starts.finditer(s, overlapped=True):
            # Must check the start is in frame
            start = match.start()
            if start % 3 == 0:
                n = s[start:]
                assert len(n) % 3 == 0, "%s is len %i" % (n, len(n))
                if strict:
                    t = translate(n, self.table)
                else:
                    # Use when missing stop codon,
                    t = "M" + translate(n[3:], self.table, to_stop=True)
                yield start, n, t  # Edited by CPT to be a generator

    def break_up_frame(self, s):
        """Returns offset, nuc, protein."""
        start = 0
        for match in self.re_stops.finditer(s, overlapped=True):
            index = match.start() + 3
            if index % 3 != 0:
                continue
            n = s[start:index]
            for (offset, n, t) in self.start_chop_and_trans(n):
                if n and len(t) >= self.min_len:
                    yield start + offset, n, t
            start = index

    def putative_genes_in_sequence(self, nuc_seq):
        """Returns start, end, strand, nucleotides, protein.
        Co-ordinates are Python style zero-based.
        """
        nuc_seq = nuc_seq.upper()
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
                answer.append((start, start - len(n), -1, n, t))
        answer.sort()
        return answer

    def get_all_peptides(self, nuc_seq):
        """Returns start, end, strand, nucleotides, protein.

        Co-ordinates are Python style zero-based.
        """
        # Refactored into generator by CPT
        full_len = len(nuc_seq)
        if self.strand != "reverse":
            for frame in range(0, 3):
                for offset, n, t in self.break_up_frame(nuc_seq[frame:]):
                    start = frame + offset  # zero based
                    yield (start, start + len(n), +1, n, t)
        if self.strand != "forward":
            rc = reverse_complement(nuc_seq)
            for frame in range(0, 3):
                for offset, n, t in self.break_up_frame(rc[frame:]):
                    start = full_len - frame - offset  # zero based
                    yield (start - len(n), start, -1, n, t)


class MGAFinder(object):

    def __init__(self, table, ftype, ends, min_len):
        self.table = table
        self.table_obj = CodonTable.ambiguous_generic_by_id[table]
        self.ends = ends
        self.ftype = ftype
        self.min_len = min_len
        self.starts = sorted(self.table_obj.start_codons)
        self.stops = sorted(self.table_obj.stop_codons)
        self.re_starts = re.compile("|".join(self.starts))
        self.re_stops = re.compile("|".join(self.stops))

    def locate(self, fasta_file, out_nuc, out_prot, out_bed, out_gff3):
        seq_format = "fasta"
        log.debug("Genetic code table %i" % self.table)
        log.debug("Minimum length %i aa" % self.min_len)

        out_count = 0

        out_gff3.write('##gff-version 3\n')

        for idx, record in enumerate(SeqIO.parse(fasta_file, seq_format)):
            for i, (f_start, f_end, f_strand, n, t) in enumerate(self.get_all_peptides(str(record.seq).upper())):
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
                    nice_strand, 0, 'ID=%s.%s.%s' % (self.ftype, idx, i + 1)])) + '\n')
        log.info("Found %i %ss", out_count, self.ftype)

    def start_chop_and_trans(self, s, strict=True):
        """Returns offset, trimmed nuc, protein."""
        if strict:
            assert s[-3:] in self.stops, s
        assert len(s) % 3 == 0
        for match in self.re_starts.finditer(s, overlapped=True):
            # Must check the start is in frame
            start = match.start()
            if start % 3 == 0:
                n = s[start:]
                assert len(n) % 3 == 0, "%s is len %i" % (n, len(n))
                if strict:
                    t = translate(n, self.table)
                else:
                    # Use when missing stop codon,
                    t = "M" + translate(n[3:], self.table, to_stop=True)
                yield start, n, t

    def break_up_frame(self, s):
        """Returns offset, nuc, protein."""
        start = 0
        for match in self.re_stops.finditer(s, overlapped=True):
            index = match.start() + 3
            if index % 3 != 0:
                continue
            n = s[start:index]
            for (offset, n, t) in self.start_chop_and_trans(n):
                if n and len(t) >= self.min_len:
                    yield start + offset, n, t
            start = index

    def putative_genes_in_sequence(self, nuc_seq):
        """Returns start, end, strand, nucleotides, protein.
        Co-ordinates are Python style zero-based.
        """
        nuc_seq = nuc_seq.upper()
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
                answer.append((start, start - len(n), -1, n, t))
        answer.sort()
        return answer

    def get_all_peptides(self, nuc_seq):
        """Returns start, end, strand, nucleotides, protein.

        Co-ordinates are Python style zero-based.
        """
        # Refactored into generator by CPT

        full_len = len(nuc_seq)
        for frame in range(0, 3):
            for offset, n, t in self.break_up_frame(nuc_seq[frame:]):
                start = frame + offset  # zero based
                yield (start, start + len(n), +1, n, t)
        rc = reverse_complement(nuc_seq)
        for frame in range(0, 3):
            for offset, n, t in self.break_up_frame(rc[frame:]):
                start = full_len - frame - offset  # zero based
                yield (start - len(n), start, -1, n, t)
