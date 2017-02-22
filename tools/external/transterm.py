#!/usr/bin/env python
import os
import sys
import argparse
import subprocess
import tempfile
import re
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF
from gff3 import feature_lambda, feature_test_type
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

REGEX_TERM = re.compile('  TERM \d+ \s* (\d+) - (\d+)\s*(\+|-) [^ ]* \s* (\d+)\s* ([0-9.-]+)\s* -([0-9.]+)\s*\|\s*(.*)')
PARTS_TERM = re.compile('  ([^ ]*)\s+([^ ]*)\s+([^ ]*)\s+([^ ]*)\s+([^ ]*)')
COLS = ("start", "end", "strand", "confidence", "hp score", "tail score",
        "notes", "5' tail", "5' stem", "loop", "3' stem", "3' tail")


def build_expterm():
    pass


def generate_annotation_file(gff3):
    # TODO: cleanup
    t = tempfile.NamedTemporaryFile(delete=False, suffix='.coords')
    for rec in GFF.parse(gff3):
        features = feature_lambda(rec.features, feature_test_type, {'type': 'CDS'}, subfeatures=False)
        for feature in sorted(features, key=lambda x: x.location.start):
            t.write('\t'.join(map(str, [
                feature.id,
                feature.location.start + 1,
                feature.location.end,
                rec.id
            ])) + '\n')
    name = t.name
    t.close()
    return name


def run_transterm(expterm, fasta, annotations):
    output = subprocess.check_output([
        'transterm', '-p', expterm, '--all-context', fasta, annotations
    ])
    return output


def pairwise(it):
    it = iter(it)
    while True:
        yield next(it), next(it)


def parse_transterm(data):
    data = data.split('SEQUENCE')[1:]
    for datum in data:
        lines = datum.split('\n')

        seq_name = lines[0][1:]
        seq_name = seq_name[0:seq_name.index(' ')]

        record = SeqRecord(
            Seq("ACTG"),
            id=seq_name
        )

        important_lines = [
            x for x in lines
            if len(x.strip()) > 0 and
            x.startswith('  ')
        ]
        # Must have an even #
        assert len(important_lines) % 2 == 0
        for idx, (a, b) in enumerate(pairwise(important_lines)):
            parsed_data = \
                REGEX_TERM.match(a).groups() + \
                PARTS_TERM.match(b).groups()

            pd = {k: v for (k, v) in zip(COLS, parsed_data)}
            # start , end  , strand , confidence , hp score , tail score , notes
            # , 5' tail         , 5'stem   , loop , 3' stem  , 3'loop
            start = int(pd['start'])
            end = int(pd['end'])

            notes = {k: pd[k] for k in COLS[6:]}
            notes2 = {}
            # However, if on - strand, we need to swap 5', 3':
            if pd['strand'] == '-':
                # Also need to flip the data
                notes2 = {
                    "[1] 5' stem": str(Seq(notes["3' stem"]).reverse_complement()),
                    "[2] loop": str(Seq(notes["loop"]).reverse_complement()),
                    "[3] 3' stem": str(Seq(notes["5' stem"]).reverse_complement()),
                    "[4] 3' tail": str(Seq(notes["5' tail"]).reverse_complement()),
                }
            else:
                notes2 = {
                    "[1] 5' stem": notes["5' stem"],
                    "[2] loop": notes["loop"],
                    "[3] 3' stem": notes["3' stem"],
                    "[4] 3' tail": notes["3' tail"],
                }

            qualifiers = {
                'score': pd['confidence'],
                'source': 'TranstermHP',
                'ID': ['terminator_%s' % idx],
            }
            current_start = min(start, end) - 1
            current_end = max(start, end)

            if pd['strand'] == '+':
                # Let's extend the current_end to include any Ts we find.
                # Take the 3' tail, and check to see how many Ts we can strip:
                #
                # Updated algo: take as many chars as possible until >1 is non-T
                # If the non-T is last, strip.
                # Otherwise (internal), leave it.
                addition = ""
                prime3tail = notes["3' tail"]
                for idx in range(len(prime3tail)):
                    addition += prime3tail[idx]
                    if addition.count('A') + addition.count('C') + addition.count('G') > 1:
                        break

                if addition[-1] != 'T':
                    addition = addition[0:-1]

                current_end += len(addition)
            else:
                addition = ""
                prime5tail = notes["5' tail"][::-1]
                for idx in range(len(prime5tail)):
                    addition += prime5tail[idx]
                    if addition.count('T') + addition.count('C') + addition.count('G') > 1:
                        break

                if addition[-1] != 'A':
                    addition = addition[0:-1]

                current_start -= len(addition)


            qualifiers.update(notes2)
            feature = SeqFeature(
                FeatureLocation(current_start, current_end),
                type="terminator",
                strand=1 if pd['strand'] == '+' else -1,
                qualifiers=qualifiers,
            )
            record.features.append(feature)

        yield record


def main(fasta, gff3, existing_expterm='', **kwargs):
    coords_file = generate_annotation_file(gff3)
    transterm_output = run_transterm(
        existing_expterm,
        fasta,
        coords_file
    )
    try:
        os.unlink(coords_file)
    except Exception:
        pass
        # Not my job

    for record in parse_transterm(transterm_output):
        yield record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export corresponding sequence in genome from GFF3', epilog="")
    parser.add_argument('fasta', help='Fasta Genome')
    parser.add_argument('gff3', help='GFF3 File')

    parser.add_argument('--min_conf', type=int, default=76, help='Only output terminators with confidence >= n')

    # parser.add_argument('--gc', type=float, default=-2.3, help='Score of a G-C pair')
    # parser.add_argument('--au', type=float, default=-0.9, help='Score of an A-U pair')
    # parser.add_argument('--gu', type=float, default=1.3, help='Score of a G-U pair')
    # parser.add_argument('--mm', type=float, default=3.5, help='Score of any other pair')
    # parser.add_argument('--gap', type=int, default=6, help='Score of a gap in the hairpin')
    # parser.add_argument('--max_hp_score', type=float, default=-2, help='Maximum allowable hairpin score')
    # parser.add_argument('--max_tail_score', type=float, default=-2.5, help='Maximum allowable tail score')
    # parser.add_argument('--max_len', type=int, default=59, help='Total extent of hairpin <= n NT long')
    # parser.add_argument('--min_stem', type=int, default=4, help='Stem must be n nucleotides long')
    # parser.add_argument('--max_loop', type=int, default=13, help='The loop portion can be no longer than n')
    # parser.add_argument('--min_loop', type=int, default=3, help='Loop portion of the hairpin must be at least n long')
    # parser.add_argument('--uwin_require', type=int, default=3, help='Number of "U" nucleotides in the --uwin_length long region.')
    # parser.add_argument('--loop_penalty', default='1,2,3,4,5,6,7,8,9,10,11', help='The cost of loops of various lengths can be set using --loop_penalty=f1,f2,f3,f4,f5,...fn, where f1 is the cost of a loop of length --min_loop, f2 is the cost of a loop of length --min_loop+1, as so on. If there are too few terms to cover up to max_loop, the last term is repeated.',)

    args = parser.parse_args()

    for record in main(existing_expterm=os.path.join(SCRIPT_PATH, 'expterm.dat'), **vars(args)):
        GFF.write([record], sys.stdout)
