#!/usr/bin/env python
import re
import argparse
import sys
import random
import itertools
import bisect
import yaml
from collections import Sequence
from Bio.Seq import Seq
from Bio import SeqIO
from itertools import groupby
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

MAXTRIES = 10000


class Enzyme(object):
    # Ugly code re-use but we don't package things nicely enough to avoid it.

    DNA_REGEX_TRANSLATIONS = {
        'A': 'A',
        'T': 'T',
        'C': 'C',
        'G': 'G',
        'N': '.',
        'M': '[AC]',
        'R': '[AG]',
        'W': '[AT]',
        'Y': '[CT]',
        'S': '[CG]',
        'K': '[GT]',

        'H': '[^G]',
        'B': '[^A]',
        'V': '[^T]',
        'D': '[^C]',
    }

    def __init__(self, forward, reverse):
        self.forward = forward
        self.reverse = reverse

    def get_regex(self):
        """Generate the regular expression objects for forward+backward cuts"""
        reg_f = self.iupac_to_regex(self.forward)
        reg_r = self.iupac_to_regex(self.reverse)
        rec_seq_f = re.compile(reg_f, re.IGNORECASE)
        rec_seq_r = re.compile(reg_r, re.IGNORECASE)
        return [rec_seq_f, rec_seq_r]

    def iupac_to_regex(self, recognition_sequence):
        """Replace IUPAC extended DNA alphabet characters with appropriate
        regular experssions from Enzyme.DNA_REGEX_TRANSLATIONS"""
        return ''.join([self.DNA_REGEX_TRANSLATIONS[x] for x in
                        recognition_sequence])


class Mutator(object):

    def __init__(self, target, mask, table=11, codondb=None, seed=42,
                 avoidCustom=None, avoidEnzyme=None, rebase=None, **kwargs):

        if seed > 0:
            random.seed(seed)

        self.masked_regions = self.parse_mask_files(mask)
        self.codon_table, self.translation_table = self.gen_opt_table(table=table)

        rebaseTmp = yaml.load(rebase)
        self.avoidance = avoidCustom if avoidCustom is not None else []

        if avoidEnzyme is None:
            avoidEnzyme = []

        for enzyme in rebaseTmp:
            if enzyme in avoidEnzyme:
                self.avoidance += Enzyme(*rebaseTmp[enzyme]).get_regex()

        # Load target data from codondb
        header = None
        lastline = False
        real_target = ':%s:' % target
        for line in codondb:
            if real_target in line:
                header = "" + line
                lastline = True
                continue

            if lastline:
                break

        self.target_table = self.parse_codondb(
            header,
            line.strip().split(' '))

    def parse_codondb(self, header, data):
        log.info("Found match: %s", header.strip())
        # Hardcoding beacuse lazy
        spsum_label = ('CGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC '
                       'UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC '
                       'GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU '
                       'CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU '
                       'AUA AUC AUU AUG UGG UAA UAG UGA')
        spsum = spsum_label.replace('U', 'T').split(' ')
        # Parse codondb/spsum
        organism_codon_usage = {k: int(v) for (k, v) in zip(spsum, data)}

        # Our real output
        target_codon_table = {}
        for amino_acid in self.codon_table:
            # Iterate over all the amino acids/codons
            possible_codons = self.codon_table[amino_acid]
            target_codon_table[amino_acid] = {}

            for codon in possible_codons:
                # Copy codon usage counts into our TCT
                target_codon_table[amino_acid][codon] = organism_codon_usage[codon]

        return target_codon_table

    def mutate(self, sequence):
        final_seq = Seq('')
        regions = self.generate_evaluatable_regions(sequence)
        if len(regions) == 1:
            final_seq += self._safeMutate(sequence)
        else:
            for (region_start, region_end, masked) in regions:
                region = sequence[region_start:region_end]
                log.info('[%s..%s %s] %s', region_start, region_end, masked, str(region.seq)[0:60])
                if masked:
                    final_seq += region
                else:
                    final_seq += self._safeMutate(region)

        return final_seq

    def _safeMutate(self, region):
        mutated_region = self._mutate(region)
        runs = 0
        while self._badRegion(mutated_region) and runs < MAXTRIES:
            mutated_region = self._mutate(region)
            runs += 1
            if runs % 100 == 0:
                log.info("... [%s] attempts", runs)
            # log.info("... [%s] %s", runs, mutated_region.seq)

        if runs >= MAXTRIES:
            log.error("Tried %s different variations, failed to find one without target sequence.", MAXTRIES)

        return mutated_region

    def _badRegion(self, sequence):
        result = False
        us = str(sequence.seq.upper())
        for query in self.avoidance:
            if isinstance(query, str):
                result = query.upper() in us
            else:
                result = len(query.findall(us)) > 0

            if result:
                return True
        return False

    def _mutate(self, seq):
        fixed_seq = ''
        for i in range(0, len(seq), 3):
            codon = str(seq[i:i + 3].seq.upper())
            if len(codon) != 3:
                fixed_seq += codon
                log.warn("Feature was not of length % 3")
                continue

            codon_aa = self.translation_table[codon]
            possible_alternates = {k: v for (k, v) in self.target_table[codon_aa].iteritems() if v > 0.1}
            if len(possible_alternates) == 0:
                raise Exception(("This should NEVER occur. It should be mathematically impossible unless "
                                 "your translation table specifies 11 possible translations for a single amino acid."))
            replacement = self.weighted_sample(possible_alternates)[0]
            log.debug('Codon: %s translates to %s, target table offers following frequencies: %s and we chose %s',
                      codon, codon_aa, possible_alternates, replacement)
            fixed_seq += replacement
        seq.seq = fixed_seq
        return seq

    def weighted_sample(self, popweights, k=1):
        return random.sample(WeightedPopulation(popweights), k=k)

    def parse_mask_files(self, bedfiles):
        if bedfiles is None:
            return []

        regions = {}
        for x in bedfiles:
            for line in x:
                bedline = line.strip().split('\t')
                chrId = bedline[0]
                if chrId not in regions:
                    regions[chrId] = []

                regions[chrId].append(tuple(map(int, bedline[1:3])))

        return regions

    def gen_opt_table(self, table=11):
        data = {}
        tntable = {}
        for codons in itertools.product('ACTG', repeat=3):
            key = ''.join(codons)
            seq = Seq(key)
            res = str(seq.translate(table=table))
            if res in data:
                data[res].append(key)
            else:
                data[res] = [key]
            tntable[key] = res
        return data, tntable

    def generate_evaluatable_regions(self, sequence):
        # Set of booleans representing masked/unmaksed regions
        regions = [True for x in range(len(sequence))]
        # for each masked in region
        if sequence.id in self.masked_regions:
            for mask in self.masked_regions[sequence.id]:
                # Loop over all values in that mask
                for i in range(mask[0], mask[1]):
                    # Set those regions to false.
                    regions[i] = False

        # now we do some post processing
        realRegions = []
        idx = 0
        # Grouping the regions by runs of Trues/Falses
        for (i, j) in groupby(regions):
            # i = the group's value, j = the group itself
            q = len(list(j))
            # Taking the length of the list, we append (start, end, masked) to
            # our realRegions
            realRegions.append((idx, idx + q, i))
            # And update our index.
            idx += q
        return realRegions


class WeightedPopulation(Sequence):
    """
    http://stackoverflow.com/questions/13047806/weighted-random-sample-in-python
    """
    def __init__(self, population_weights):
        self.population = population_weights.keys()
        self.cumweights = []
        cumsum = 0  # compute cumulative weight
        for weight in population_weights.itervalues():
            cumsum += weight
            self.cumweights.append(cumsum)

    def __len__(self):
        return self.cumweights[-1]

    def __getitem__(self, i):
        if not 0 <= i < len(self):
            raise IndexError(i)
        return self.population[bisect.bisect(self.cumweights, i)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('sequence', type=file, help='sequence file')
    parser.add_argument('target', type=str, help='target organism')
    parser.add_argument('--table', type=int, help='Translation table #', default=1)
    parser.add_argument('--mask', type=file, nargs='*', help='Regions to mask for mutations')
    parser.add_argument('--codondb', type=file, help='Average codon database')
    parser.add_argument('--rebase', type=file, help='Rebase yaml file')
    parser.add_argument('--seed', type=int, help='Random seed. 0 means choose randomly at runtime', default=0)
    parser.add_argument('--avoidCustom', nargs='*', help='Sequences to avoid')
    parser.add_argument('--avoidEnzyme', nargs='*', help='Enzymes to avoid')

    args = parser.parse_args()
    m = Mutator(**vars(args))

    for sequence in SeqIO.parse(args.sequence, 'fasta'):
        SeqIO.write(m.mutate(sequence), sys.stdout, 'fasta')
