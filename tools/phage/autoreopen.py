#!/usr/bin/env argparse
import os
import json
import argparse
import subprocess
from enum import Enum
from Bio import SeqIO
from BCBio import GFF
from cpt_convert_mga_to_gff3 import mga_to_gff3
from gff3 import feature_lambda
from safe_reopen import extract_gff3_regions, gaps
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, 'reopen-data')
# TODO: This tool depends on PY2K ONLY TOOLS. THIS WILL(MAY?) NOT FUNCTINO UNDER PY3.
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('autoreopen')


class AutoreopenException(Exception):
    """Class to base exceptions off of. This allows us to catch normal program
    exceptions without catching worse things as well."""
    pass

class Evidence(Enum):
    NONE = 0
    Assumption = 1
    BLAST = 2
    PhageTerm = 3


class PhageType(Enum):
    ShortTR= 1
    LongTR = 2
    Prime3Cos = 3
    Prime5Cos = 4
    PacHeadful = 6
    T4Headful = 7
    MuLike = 8
    Unknown = 9


class PhageReopener:

    def __init__(self, fasta, fastq1, fastq2, closed=False):
        self.rec_file = fasta
        self.fasta = SeqIO.read(fasta, 'fasta')
        self.closed = closed
        self.fq1 = fastq1
        self.fq2 = fastq2

    def _orfCalls(self):
        fnmga = os.path.join(DATA_DIR, str(self.fasta.id) + '.mga')
        if not os.path.exists(fnmga):
            # Run MGA
            subprocess.check_call([
                'mga_linux_x64', '-s', self.rec_file.name, '>', fnmga
            ], shell=True)

        # Convert to gff3
        fn = str(self.fasta.id) + '.mga.gff3'
        self.mga_gff3 = fn
        with open(fnmga, 'r') as handle, open(fn, 'w') as output:
            self.rec_file.seek(0)
            for result in mga_to_gff3(handle, self.rec_file):
                # Store gFF3 data in self in order to access later.
                self.mga_rec = result
                GFF.write([result], output)

        # Process a feature id -> feature table in mem.
        self.featureDict = {}
        for f in feature_lambda(self.mga_rec.features, lambda x: True, {}, subfeatures=True):
            self.featureDict[f.qualifiers['ID'][0]] = f

        # Extract
        fnfa = str(self.fasta.id) + '.mga.fa'
        subprocess.check_call([
            'python2', os.path.join(os.pardir, 'gff3', 'gff3_extract_sequence.py'),
            '--feature_filter', 'CDS',
            self.rec_file.name, fn,
        ], stdout=open(fnfa, 'w'))

        # Translate
        fnpfa = str(self.fasta.id) + '.mga.pfa'
        subprocess.check_call([
            'python2', os.path.join(os.pardir, 'fasta', 'fasta_translate.py'),
            '--table', '11', '--strip_stops', '--target', 'protein',
            fnfa
        ], stdout=open(fnpfa, 'w'))
        return fnpfa

    def _runPhageTerm(self):
        fn = os.path.join(DATA_DIR, str(self.fasta.id) + '.json')
        if not os.path.exists(fn):
            subprocess.check_call([
                'python2', os.path.join(os.pardir, 'external', 'phageterm', 'PhageTerm.py'),
                '-c', '6',
                '-f', self.fq1.name,
                '-n', self.fasta.id,
                '-r', self.rec_file.name,
                '-p', self.fq2.name,
                '-s', '30'
            ])

        with open(fn, 'r') as handle:
            return json.load(handle)

    def _analysePhageTerm(self, results):
        phtype = None
        opening = results['reopening']

        if results['P_class'] == "COS (3')":
            phtype = PhageType.Prime3Cos
        elif results['P_class'] == "COS (5')":
            phtype = PhageType.Prime5Cos
        elif results['P_class'] == 'DTR (long)':
            phtype = PhageType.LongTR
        elif results['P_class'] == 'DTR (short)':
            phtype = PhageType.ShortTR
        elif results['P_class'] == 'Headful (pac)':
            phtype = PhageType.PacHeadful
        elif results['P_class'] == 'Mu-like':
            phtype = PhageType.MuLike
        else:
            raise AutoreopenException("Cannot interpret this automatically, yet.")

        return (phtype, opening, Evidence.PhageTerm)

    def _safeOpeningLocationForFeature(self, feature):
        """Given a feature, find a 'safe' location to re-open the genome.

        Safe here means not in the middle of a gene. Thus, we search upstream
        in the genome until we find a place not covered by genes.
        """
        occupied_regions = extract_gff3_regions([self.mga_gff3])
        # Only acceptable places to open a genome.
        gaps_in_data = gaps(occupied_regions[self.fasta.id])

        best = None
        for gap in gaps_in_data:
            if feature.location.strand > 0:
                # If the end of the gap is upstream
                if gap[1] < feature.location.start:
                    best = gap
            else:
                # if beginning of gap is upstream
                if gap[0] > feature.location.end:
                    best = gap
        # Ok, we /should/ have a 'best' gap (But I think I'm missing an edge
        # case since this isn't a circular data structure...)
        # Next, get midpoint
        mid = sum(best) / 2
        return mid

    def _blast(self, pfa):
        blast_results = subprocess.check_output([
            'blastp', '-query', pfa,
            '-db', os.path.join(SCRIPT_DIR, 'test-data', 'merged'),
            '-outfmt', '6 sseqid evalue pident qseqid',
            '-evalue', '0.001',
        ]).split('\n')

        def blast_filter(d):
            if len(d.strip()) == 0:
                return False
            q = d.split('\t')
            q[2] = float(q[2])
            return q[2] > 50

        # Filter only good quality hits (>50%)
        blast_results = filter(blast_filter, blast_results)
        # Pull out the type
        blast_result = []
        for hit in blast_results:
            (sseqid, evalue, pident, qseqid) = hit.split('\t')
            blast_result.append((sseqid[:sseqid.index('_')], qseqid))

        # reduce to set
        blast_result = list(set(blast_result))
        ph_type = None
        # Handle results
        if len(blast_result) == 1:
            blast_type, protein_name = blast_result[0]

            if blast_type == '3primeCos':
                ph_type = PhageType.Prime3Cos
            elif blast_type == '5primeCos':
                ph_type = PhageType.Prime5Cos
            elif blast_type == 'long-tr':
                ph_type = PhageType.LongTR
            elif blast_type == 'short-tr':
                ph_type = PhageType.ShortTR
            elif blast_type == 'pac-headful':
                ph_type = PhageType.PacHeadful
            elif blast_type == 't4-headful':
                ph_type = PhageType.T4Headful
            else:
                log.warning("PhageType %s??", blast_type)
                return None

            f = self.featureDict[protein_name]
            return ph_type, ('Forward' if f.strand > 0 else 'Reverse', self._safeOpeningLocationForFeature(f))
        else:
            log.warning("%s blast results for this protein", len(blast_results))
            return None

    def detectType(self):
        # Naive annotation, we run these NO MATTER WHAT.
        self.protein_fasta = self._orfCalls()
        (blast_type, blast_reopen_location) = self._blast(self.protein_fasta)

        # Try their tool
        try:
            results = self._runPhageTerm()
            (phageterm_type, phageterm_location) = self._analysePhageTerm(results)
            return (phageterm_type, phageterm_location, Evidence.PhageType)
        except AutoreopenException:
            # Here we've failed to get a useful result from Phage Term, so
            # continue on.
            pass

        # Next we failover to blast results.
        if blast_type:
            # This will give us a tuple, (PhageType, locationToReopen)
            return (blast_type, blast_reopen_location, Evidence.BLAST)

        # Failing that, if the genome is closed
        if self.closed:
            # we assume it is headful
            return (PhageType.PacHeadful, None, Evidence.Assumption)
        else:
            # Otherwise we truly have no clue
            return (PhageType.Unknown, None, Evidence.NONE)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Automatic re-opening tool')
    parser.add_argument('fasta', type=argparse.FileType("r"))
    parser.add_argument('fastq1', type=argparse.FileType("r"))
    parser.add_argument('fastq2', type=argparse.FileType("r"))
    parser.add_argument('--closed', action='store_true')
    args = parser.parse_args()

    # create new IntronFinder object based on user input
    pr = PhageReopener(**vars(args))
    print(pr.detectType())
