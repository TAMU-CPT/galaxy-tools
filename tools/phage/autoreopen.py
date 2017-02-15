#!/usr/bin/env python
import os
import re
import json
import shutil
import argparse
import subprocess
import numpy as np
from enum import Enum
from Bio import SeqIO
from BCBio import GFF
from cpt_convert_mga_to_gff3 import mga_to_gff3
from gff3 import feature_lambda
from safe_reopen import extract_gff3_regions, gaps
from jinja2 import Environment, FileSystemLoader
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, 'reopen-data')
# TODO: This tool depends on PY2K ONLY TOOLS. THIS WILL(MAY?) NOT FUNCTINO UNDER PY3.
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('autoreopen')
env = Environment(loader=FileSystemLoader(SCRIPT_DIR), trim_blocks=True, lstrip_blocks=True)


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


class BlastBasedReopening(object):
    """Re-open genomes based on blast"""

    def __init__(self, db='test-data/canonical_nucl', genome=None, genome_real=None, protein=False):
        self.genome = genome
        self.genome_nucleotide_real = genome_real
        self.protein = protein
        self.db = db # 'test-data/canonical_nucl'
        self.canonical_map = {
        'NC_025442.1': 'Chi',
        'NC_005880.2': 'Phage K',
        'NC_008720.1': 'N4',
        'NC_004775.1': 'epsilon15',
        'NC_005282.1': 'Felix01',
        'NC_001416.1': 'lambda',
        'NC_000929.1': 'Mu',
        'NC_001901.1': 'N15',
        'NC_005856.1': 'P1',
        'NC_002371.2': 'P22',
        'NC_001895.1': 'P2',
        'NC_001609.1': 'P4',
        'NC_011048.1': 'phi29',
        'NC_005045.1': 'phiKMV',
        'NC_004831.2': 'SP6',
        'NC_011421.1': 'SPO1',
        'NC_004166.2': 'SPP1',
        'NC_005833.1': 'T1',
        'NC_003298.1': 'T3',
        'NC_000866.4': 'T4',
        'NC_005859.1': 'T5',
        'NC_001604.1': 'T7',
        }

    def blast2Xmfa(self, blast_results, name='nucl'):
        tpl = env.get_template('plot_template.html')
        with open('%s.html' % name, 'w') as html:
            html.write(tpl.render(filename='%s.json' % name).encode('utf-8'))

        with open('%s.json' % name, 'w') as handle:
            json.dump({
                "xmfa": "%s_regions.json" % name,
                "fasta": [
                    {
                        "path": "n.a",
                        "name": "SubjectGenome",
                        "length": max([max(x['sstart'], x['send']) for x in blast_results]),
                    },
                    {
                        "path": "n.a",
                        "name": "QueryGenome",
                        "length": max([max(x['qend'], x['qstart']) for x in blast_results]),
                    }
                ],
                "gff3": ["n.a", "n.a"]
            }, handle, sort_keys=True, indent=2)

        with open('%s_regions.json' % name, 'w') as handle:
            region_data = []
            for hit in blast_results:
                region_data.append([
                    {
                        "rid" : "1:fake",
                        "comment" : "",
                        "start" : hit['qstart'],
                        "seq" : "",
                        "id" : "1",
                        "strand" : 1 if hit['sstart'] > hit['send'] else -1,
                        "end" : hit['send']
                    },
                    {
                        "rid" : "2:fake",
                        "comment" : "",
                        "start" : hit['qstart'],
                        "seq" : "",
                        "id" : "2",
                        "strand" : 1 if hit['qstart'] > hit['qend'] else -1,
                        "end" : hit['send']
                    },
                ])

            json.dump(region_data, handle, sort_keys=True, indent=2)

    def evaluate(self):
        if self.protein:
            blast_results = self.getBlastP(self.genome)
            self.blast2Xmfa(blast_results, 'prot')
        else:
            blast_results = self.getBlastN(self.genome)
            self.blast2Xmfa(blast_results, 'nucl')


        orientation = '+'
        current_genome_file = None
        if not self.protein:
            # If not a protein file, we check the orientation. We probably
            # should do this for protein too, but ... yikes.
            genome_seq = SeqIO.read(self.genome, 'fasta')
            current_genome_file = self.genome

            # If more than half the hits are 'backwards'
            if self.shouldReverse(blast_results) > 0.5:
                orientation = '-'
                new_genome_file = self.genome.replace('.fa', '.revcom.fa')
                rec = genome_seq.reverse_complement(id=True, description=True)
                rec.description += ' autoreopen.revcom'
                SeqIO.write([rec],
                            new_genome_file, 'fasta')
                current_genome_file = new_genome_file
                blast_results = self.getBlastN(new_genome_file)

        # TODO: Filter out just ONE genome + best results
        # blast_results = blast_results

        # ref_genome_hits = [self.avg(hit['sstart'], hit['send']) for hit in blast_results]
        # our_genome_hits = [self.avg(hit['qstart'], hit['qend']) for hit in blast_results]
        # # Now we sort relative to ref genome
        # hits = zip(our_genome_hits, ref_genome_hits)
        # hits = sorted(hits, key=lambda x: x[1])
        # # Then we re-construct *_genome_hits from the sorted/zipped version
        # ref_genome_hits = [x[1] for x in hits]
        # our_genome_hits = [x[0] for x in hits]
        # These are now two independent lists, sorted relative to ref genome. We
        # can shift our_genome_hits without affecting ref genome ordering.
        # (score, (index, m, c)) = self.scoreResults(ref_genome_hits, our_genome_hits)

        # Now on to plotting the real data.
        totalData = [(hit['qstart'], hit['qend'], hit['sstart'], hit['send']) for hit in blast_results]
        totalData = sorted(totalData, key=lambda x: self.avg(x[2], x[3]))

        if self.protein:
            # In order to use this we need to re-open according to our new end,
            # and then re-call genes / export to pfa / re-blast. This is SUPER
            # ugly and SUPER hacky.
            new_genome_file_reopened = self.genome_nucleotide_real
            customPhageReopener = PhageReopener(open(self.genome_nucleotide_real, 'r'), None, None)
            protein_fasta = customPhageReopener._orfCalls()
            # This gets a new protein_fasta file, which we can re-blast and re-analyse
            updated_blast_results = self.getBlastP(protein_fasta)
            self.blast2Xmfa(updated_blast_results, 'prot2')
        else:
            genome_seq = SeqIO.read(current_genome_file, 'fasta')
            new_genome_file_reopened = self.genome.replace('.fa', '.reopened.fa')
            with open(new_genome_file_reopened, 'w') as handle:
                loc = totalData[0][0]
                rec = genome_seq[loc:] + genome_seq[0:loc]
                rec.description += ' autoreopen.reopen(%s)' % loc
                SeqIO.write([rec], handle, 'fasta')

            updated_blast_results = self.getBlastN(new_genome_file_reopened)
            self.blast2Xmfa(updated_blast_results, 'nucl2')


        return {
            'location': totalData[0],
            'orientation': orientation,
            'results': blast_results,
            'reopened_file': new_genome_file_reopened,
        }

    def getBlastN(self, path):
        blast_results = subprocess.check_output([
            'blastn', '-db', self.db,
            '-query', path, '-outfmt', '6 sseqid evalue pident qstart qend sstart send'
        ]).strip().split('\n')
        blast_results = [x.split('\t') for x in blast_results]
        blast_results = [
            {
                'sseqid': self.canonical_map[sseqid],
                'evalue': float(evalue),
                'pident': float(pident),
                'qstart': int(qstart),
                'qend': int(qend),
                'sstart': int(sstart),
                'send': int(send),
            }
            for (sseqid, evalue, pident, qstart, qend, sstart, send) in blast_results
            if float(evalue) < 0.001 and float(pident) > 80
        ]

        sseqids = [x['sseqid'] for x in blast_results]
        if len(set(sseqids)) > 1:
            log.warning("DO not know how to handle. This WILL MISBEHAVE")
        return blast_results

    def getBlastP(self, path):
        print('getBlastP', path)
        cmd = [
            'blastp', '-db', self.db,
            '-query', path, '-outfmt', '6 qseqid sseqid evalue pident qstart qend sstart send'
        ]
        print(' '.join(cmd))
        blast_results = subprocess.check_output(cmd).strip().split('\n')
        blastp_locmap = open(os.path.join(SCRIPT_DIR, 'test-data', 'merged.pfa.idmap'), 'r').read().strip().split('\n')

        def extremes(data):
            mi = None
            ma = None

            for d in data:
                d = int(d)
                if mi is None:
                    mi = d
                if ma is None:
                    ma = d

                ma = max(ma, d)
                mi = min(mi, d)
            return mi, ma

        blastp_locmap2 = {}
        for x in blastp_locmap:
            (key, value) = x.split('\t')
            (mi, ma) = extremes([y for y in re.split('[^0-9]+', value) if len(y) > 0])
            #print([y for y in re.split('[^0-9]+', value) if len(y) > 0])
            if 'complement' in value:
                blastp_locmap2[key] = (ma, mi)
            else:
                blastp_locmap2[key] = (mi, ma)

        #import sys; sys.exit()
        input_locmap = {}
        for tmprec in SeqIO.parse(path, 'fasta'):
            desc = tmprec.description
            locstr = desc[desc.index('Location=') + 10:-2]
            locdata = re.split('[^0-9+-]+', locstr)

            if locdata[2] == '+':
                input_locmap[tmprec.id] = (int(locdata[0]), int(locdata[1]))
            else:
                input_locmap[tmprec.id] = (int(locdata[1]), int(locdata[0]))

        blast_results = [x.split('\t') for x in blast_results]
        blast_results = [
            {
                'qseqid': qseqid[qseqid.rindex('_') + 1:],
                'sseqid': self.canonical_map[sseqid[4:sseqid.index('_prot_')]],
                'evalue': float(evalue),
                'pident': float(pident),
                'qstart': input_locmap[qseqid][0],
                'qend': input_locmap[qseqid][1],
                'sstart': blastp_locmap2[sseqid][0],
                'send': blastp_locmap2[sseqid][1],
            }
            for (qseqid, sseqid, evalue, pident, qstart, qend, sstart, send) in blast_results
            if float(evalue) < 0.001 and float(pident) > 70
        ]

        sseqid_set = set([x['sseqid'] for x in blast_results])
        if len(sseqid_set) > 1:
            log.warning("DO not know how to handle. This WILL MISBEHAVE. " + ','.join(sseqid_set))
        return blast_results

    def shouldReverse(self, blast_results):
        should_reverse = 0
        total = len(blast_results)

        for hit in blast_results:
            qdir = 1 if hit['qstart'] > hit['qend'] else -1
            sdir = 1 if hit['sstart'] > hit['send'] else -1
            if (qdir > 0 and sdir < 0) or (qdir < 0 and sdir > 0):
                should_reverse += 1

        return float(should_reverse) / total

    @classmethod
    def rotate(cls, data, index):
        return data[index:] + data[0:index]

    @classmethod
    def avg(cls, a, b):
        return int(float(a + b) / 2)

    def scoreResults(self, ref_genome_hits, our_genome_hits):

        best_score = None
        best_params = None

        # This no longer needs to be re-sorted.
        for i in range(len(our_genome_hits)):
            hits = zip(self.rotate(our_genome_hits, i), ref_genome_hits)
            x_position = np.array([hit[0] for hit in hits])
            y_position = np.array([hit[1] for hit in hits])
            (m, cx), r2, c, d, e = np.polyfit(x_position, y_position, 1, full=True)
            r2 = r2[0]

            # Now we look to minimize R**2. Could also optimize for m == 1?
            if best_score is None:
                best_score = r2
                best_params = (i, m, cx)

            if r2 < best_score:
                best_score = r2
                best_params = (i, m, cx)

        return best_score, best_params


class PhageReopener:

    def __init__(self, fasta, fastq1, fastq2, closed=False, data_dir=None):
        # fasta, fastq1, fastq2 are all file handles
        self.base_name = fasta.name.replace('.fa', '')
        self.rec_file = fasta
        self.data_dir = data_dir
        self.fasta = SeqIO.read(fasta, 'fasta')
        self.fq1 = fastq1
        self.fq2 = fastq2
        self.closed = closed

    def _orfCalls(self):
        fnmga = self.base_name + '.mga'
        print(fnmga)
        if not os.path.exists(fnmga):
            # Run MGA
            subprocess.check_call([
                'mga_linux_x64', '-s', self.rec_file.name
            ], stdout=open(fnmga, 'w'))

        # Convert to gff3
        fn = self.base_name + '.mga.gff3'
        self.mga_gff3 = fn
        with open(fnmga, 'r') as handle, open(fn, 'w') as output:
            self.rec_file.seek(0)
            print(self.rec_file.name)
            for result in mga_to_gff3(handle, self.rec_file):
                # Store gFF3 data in self in order to access later.
                self.mga_rec = result
                GFF.write([result], output)

        # Process a feature id -> feature table in mem.
        self.featureDict = {}
        for f in feature_lambda(self.mga_rec.features, lambda x: True, {}, subfeatures=True):
            self.featureDict[f.qualifiers['ID'][0]] = f

        # Extract
        fnfa = self.base_name + '.mga.fa'
        self.fnfa = fnfa
        subprocess.check_call([
            'python2', os.path.join(os.pardir, 'gff3', 'gff3_extract_sequence.py'),
            '--feature_filter', 'CDS',
            self.rec_file.name, fn,
        ], stdout=open(fnfa, 'w'))

        # Translate
        fnpfa = self.base_name +  '.mga.pfa'
        self.fnpfa = fnpfa
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
            phtype = PhageType.Unknown

        if phtype == PhageType.Unknown:
            opening = ['Unknown']

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
        ph_type = PhageType.Unknown
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
                # log.warning("PhageType %s??", blast_type)
                ph_type = PhageType.Unknown

            f = self.featureDict[protein_name]
            return blast_results, ph_type, ['forward' if f.strand > 0 else 'reverse', self._safeOpeningLocationForFeature(f)]
        else:
            log.warning("%s blast results for this protein", len(blast_results))
            return blast_results, ph_type, ['Unknown']

    def giveEvidence(self):
        kwargs = {}
        # Naive annotation, we run these NO MATTER WHAT.
        self.protein_fasta = self._orfCalls()
        (blast_results, blast_type, blast_reopen_location) = self._blast(self.protein_fasta)

        # Try their tool
        results = self._runPhageTerm()
        (phageterm_type, phageterm_location, phageterm_evidence) = self._analysePhageTerm(results)
        phageterm = {'type': phageterm_type.name,
                     'location': phageterm_location}
        kwargs['PhageTerm'] = phageterm

        # Next we failover to blast results of terminases.
        blast = {'type': blast_type.name,
                 'location': blast_reopen_location,
                 'results': [x.split() for x in blast_results]}
        kwargs['BLAST'] = blast

        # Then we do alignment relative to canonical
        bbr_nucl = BlastBasedReopening(genome=self.rec_file.name, db='test-data/canonical_nucl', protein=False)
        kwargs['canonical_nucl'] = bbr_nucl.evaluate()
        kwargs['canonical_nucl']['location'] = kwargs['canonical_nucl']['location'][0]
        ## location, results, orientation

        bbr_prot = BlastBasedReopening(genome=self.fnpfa, db='test-data/canonical_prot', protein=True,
                                       genome_real=kwargs['canonical_nucl']['reopened_file'])
        kwargs['canonical_prot'] = bbr_prot.evaluate()
        # if we reversed in canonical nucl, we'll reverse the protein results
        # as well to make them match up better.
        if kwargs['canonical_nucl']['orientation'] == '-':
            kwargs['canonical_prot']['orientation'] = '-'
            kwargs['canonical_prot']['location'] = kwargs['canonical_prot']['location'][1]

        # Failing that, if the genome is closed
        if self.closed:
            # we assume it is headful
            # yield (PhageType.PacHeadful, None, Evidence.Assumption)
            pass
        else:
            # Otherwise we truly have no clue
            # yield (PhageType.Unknown, None, Evidence.NONE)
            pass

        tpl = env.get_template('template.html')
        with open('report.html', 'w') as html:
            html.write(tpl.render(**kwargs).encode('utf-8'))


def main():
    parser = argparse.ArgumentParser(description='Automatic re-opening tool')
    parser.add_argument('fasta', type=argparse.FileType("r"))
    parser.add_argument('fastq1', type=argparse.FileType("r"))
    parser.add_argument('fastq2', type=argparse.FileType("r"))
    parser.add_argument('--closed', action='store_true')
    parser.add_argument('--data_dir', help='Directory for HTML files to go.')
    args = parser.parse_args()

    # create new IntronFinder object based on user input
    pr = PhageReopener(**vars(args))
    pr.giveEvidence()


if __name__ == '__main__':
    main()
