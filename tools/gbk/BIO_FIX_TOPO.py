import Bio.GenBank


def record_end(self, content):
    """Clean up when we've finished the record.
    """
    from Bio import Alphabet
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq, UnknownSeq

    # Try and append the version number to the accession for the full id
    if not self.data.id:
        assert 'accessions' not in self.data.annotations, self.data.annotations['accessions']
        self.data.id = self.data.name  # Good fall back?
    elif self.data.id.count('.') == 0:
        try:
            self.data.id += '.%i' % self.data.annotations['sequence_version']
        except KeyError:
            pass

    # add the sequence information
    # first, determine the alphabet
    # we default to an generic alphabet if we don't have a
    # seq type or have strange sequence information.
    seq_alphabet = Alphabet.generic_alphabet

    # now set the sequence
    sequence = "".join(self._seq_data)

    if self._expected_size is not None and len(sequence) != 0 and self._expected_size != len(sequence):
        import warnings
        from Bio import BiopythonParserWarning
        warnings.warn("Expected sequence length %i, found %i (%s)." %
                      (self._expected_size, len(sequence), self.data.id),
                      BiopythonParserWarning)

    if self._seq_type:
        # mRNA is really also DNA, since it is actually cDNA
        if 'DNA' in self._seq_type.upper() or 'MRNA' in self._seq_type.upper():
            seq_alphabet = IUPAC.ambiguous_dna
        # are there ever really RNA sequences in GenBank?
        elif 'RNA' in self._seq_type.upper():
            # Even for data which was from RNA, the sequence string
            # is usually given as DNA (T not U).  Bug 2408
            if "T" in sequence and "U" not in sequence:
                seq_alphabet = IUPAC.ambiguous_dna
            else:
                seq_alphabet = IUPAC.ambiguous_rna
        elif 'PROTEIN' in self._seq_type.upper() or self._seq_type == "PRT":  # PRT is used in EMBL-bank for patents
            seq_alphabet = IUPAC.protein  # or extended protein?
        # work around ugly GenBank records which have circular or
        # linear but no indication of sequence type
        elif self._seq_type in ["circular", "linear", "unspecified"]:
            pass
        # we have a bug if we get here
        else:
            raise ValueError("Could not determine alphabet for seq_type %s" % self._seq_type)

        # Also save the chomosome layout
        if 'circular' in self._seq_type.lower():
            self.data.annotations['topology'] = 'circular'
        elif 'linear' in self._seq_type.lower():
            self.data.annotations['topology'] = 'linear'

    if not sequence and self.__expected_size:
        self.data.seq = UnknownSeq(self._expected_size, seq_alphabet)
    else:
        self.data.seq = Seq(sequence, seq_alphabet)

Bio.GenBank._FeatureConsumer.record_end = record_end
