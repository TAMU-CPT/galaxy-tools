# Copyright 2015-2015 by Eric Rasche.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.AlignIO support for "xmfa" output from Mauve/ProgressiveMauve.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).
"""

from __future__ import print_function

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from .Interfaces import AlignmentIterator, SequentialAlignmentWriter

__docformat__ = "restructuredtext en"


class MauveIterator(AlignmentIterator):
    """Mauve xmfa alignment iterator."""

    _header = None  # for caching lines between __next__ calls

    def __next__(self):
        handle = self.handle

        if self._header is None:
            line = handle.readline()
        else:
            # Header we saved from when we were parsing
            # the previous alignment.
            line = self._header
            self._header = None

        if not line:
            raise StopIteration

        # Whitelisted headers we know about
        if not line.strip().split()[0].startswith('#FormatVersion'):
            raise ValueError("XMFA files should start with #FormatVersion")
        else:
            version = line.strip().split()[1]

        # There should be two blank lines after the header line
        line = handle.readline()
        while line.strip().startswith('#'):
            line = handle.readline()

        # If the alignment contains entries with the same sequence
        # identifier (not a good idea - but seems possible), then this
        # dictionary based parser will merge their sequences.  Fix this?
        ids = []
        seqs = {}
        consensus = ""
        seq_cols = None  # Used to extract the consensus
        passed_end_alignment = False

        latest_id = None
        while True:
            line = handle.readline()
            if not line:
                break  # end of file
            line = line.strip()

            if line == '=':
                # There may be more data, but we've reached the end of this
                # alignment
                passed_end_alignment = True
            elif line.startswith('>'):
                assert not passed_end_alignment
                parts = line.split()
                id, region = parts[1].split(':')

                if id not in ids:
                    ids.append(id)

                seqs.setdefault(id, '')
                latest_id = id
            else:
                assert not passed_end_alignment
                if latest_id is None:
                    raise ValueError("Saw sequence before definition line")
                seqs[latest_id] += line

        assert len(seqs) <= len(ids)

        self.ids = ids
        self.sequences = seqs

        if ids and seqs:
            alignment_length = len(list(seqs.values())[0])
            records = []
            for id in ids:
                seq = seqs[id]
                if alignment_length != len(seq):
                    raise ValueError("Sequences have different lengths, or repeated identifier")
                name, start, end = self._identifier_split(id)
                record = SeqRecord(Seq(seq, self.alphabet), id=id, name=name, description=id)

                if start is not None:
                    record.annotations["start"] = int(start)

                if end is not None:
                    record.annotations["end"] = int(end)

                records.append(record)
            alignment = MultipleSeqAlignment(records, self.alphabet)
            return alignment
        else:
            raise StopIteration

    def _identifier_split(self, identifier):
        """Returns (name, start, end) string tuple from an identifier"""
        id, loc = identifier.split(':')
        start, end = loc.split('-')
        return id, start, end
