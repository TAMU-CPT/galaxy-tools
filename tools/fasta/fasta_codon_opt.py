#!/usr/bin/env python
import argparse
import sys
import random
import itertools
import bisect

from collections import Sequence
from Bio.Seq import Seq
from Bio import SeqIO
from operator import itemgetter
from itertools import groupby
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


class Mutator(object):

    def __init__(self, sequence, target, mask, table=11, codondb=None, seed=42):
        if seed > 0:
            random.seed(seed)

        self.sequence = SeqIO.read(sequence, 'fasta')
        self.masked_regions = self.parse_mask_files(mask, filterId=self.sequence.id)
        self.codon_table, self.translation_table = self.gen_opt_table(table=table)

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

    def mutate(self):
        final_seq = Seq('')
        regions = self.generate_evaluatable_regions()
        if len(regions) == 1:
            final_seq += self._mutate(self.sequence)
        else:
            for (region_start, region_end, masked) in regions:
                region = self.sequence[region_start:region_end]
                if masked:
                    final_seq += region
                else:
                    final_seq += self._mutate(region)

        return final_seq

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

    def parse_mask_files(self, bedfiles, filterId=None):
        if bedfiles is None:
            return []

        regions = []
        for x in bedfiles:
            for line in x:
                bedline = line.strip().split('\t')
                chrId = bedline[0]
                if chrId == filterId:
                    regions.append(tuple(map(int, bedline[1:3])))

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

    def generate_evaluatable_regions(self):
        regions = [(0, len(self.sequence), True)]
        for mask in self.masked_regions:
            # Figure out which regions are overlapped and must be split (There /should/only be 1 or 0)
            overlapping = []
            nonoverlapping = []
            for region in regions:
                if self._overlap(region, mask):
                    overlapping.append(region)
                else:
                    nonoverlapping.append(region)

            # If there's an overlapping region, we split it.
            splitregions = []
            for region in overlapping:
                intersection = self.contiguous(self._overlap(region, mask))[0]
                ri = (intersection[0], intersection[1], False)
                left = (region[0], intersection[0], True)
                right = (intersection[1], region[1], True)

                if left[1] - left[0] != 0:
                    splitregions.append(left)

                if right[1] - right[0] != 0:
                    splitregions.append(right)

                splitregions.append(ri)

            # Finally we correct our existing region set with the non-overlapping regions, and our new split regions
            regions = nonoverlapping + splitregions
        return sorted(regions, key=lambda x: x[0])

    def contiguous(self, data):
        # http://stackoverflow.com/a/2154437
        ranges = []
        for key, group in groupby(enumerate(data), lambda (index, item): index - item):
            group = map(itemgetter(1), group)
            ranges.append((group[0], group[-1] + 1))
        return ranges

    def _overlap(self, a, b):
        return set(range(*a[0:2])).intersection(range(*b[0:2]))


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
    parser.add_argument('--mask', type=file, nargs='*', help='Regions to mask from mutations')
    parser.add_argument('--codondb', type=file, help='Average codon database')
    parser.add_argument('--seed', type=int, help='Random seed. 0 means choose randomly at runtime', default=0)

    args = parser.parse_args()
    m = Mutator(**vars(args))

    SeqIO.write(m.mutate(), sys.stdout, 'fasta')
