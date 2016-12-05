#!/usr/bin/env python
from Bio.SeqUtils import GC_skew
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import argparse
import logging
logging.basicConfig()
log = logging.getLogger()


def _req_dna(f):

    def f2(*args, **kwargs):
        if kwargs['seqtype'] == 'protein':
            log.warn("Cannot calculate for proteins")
            return 0
        return f(*args, **kwargs)
    return f2


def _translate(f):

    def new_f(*args, **kwargs):
        if kwargs['seqtype'] == 'dna':
            args = args[0].seq.translate(table=kwargs.get('table', 11), cds=True),
            kwargs['seqtype'] = 'protein'
        return f(*args, **kwargs)
    return new_f


def _sequence(f):

    def str_seq(*args, **kwargs):
        args = str(args[0].seq),
        return f(*args, **kwargs)

    return str_seq


def _nostop(f):

    def remove_stop(*args, **kwargs):
        if kwargs['seqtype'] == 'protein':
            if args[0].endswith('*'):
                args = args[0][0:-1],
        else:
            args = args[0][0:-3],

        return f(*args, **kwargs)

    return remove_stop


@_translate
@_sequence
def aromaticity(sequence, **kwargs):
    return ProteinAnalysis(sequence).aromaticity()


@_sequence
@_req_dna
def gc(sequence, **kwargs):
    return GC(sequence)


@_sequence
@_req_dna
def gc_skew(sequence, **kwargs):
    values = GC_skew(sequence)
    # Average GC skew
    if len(values) == 0:
        return 0
    return reduce(lambda x, y: x + y, values) / float(len(values))


@_translate
@_sequence
def iep(sequence, **kwargs):
    return ProteinAnalysis(sequence).isoelectric_point()


@_translate
@_sequence
@_nostop
def instability(sequence, **kwargs):
    return ProteinAnalysis(sequence).instability_index()


@_translate
@_sequence
@_nostop
def mw(sequence, **kwargs):
    return ProteinAnalysis(sequence).molecular_weight()


@_translate
@_sequence
@_nostop
def length(sequence, **kwargs):
    return len(sequence)


@_req_dna
@_sequence
def start_codon(sequence, **kwargs):
    return sequence[0:3].upper()


@_req_dna
@_sequence
def stop_codon(sequence, **kwargs):
    return sequence[-3:].upper()


def main(fasta, func, protein):
    fn = globals()[func]
    print '# ID\tvalue'
    for record in SeqIO.parse(fasta, 'fasta'):
        print '{0.id}\t{1}'.format(
            record,
            fn(record, seqtype='protein' if protein else 'dna')
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sequence Properties')

    funcs = [x for x in locals().keys()
             if not x[0].upper() == x[0] and x[0] != '_' and x not in ('log', 'main', 'argparse', 'logging', 'parser')]

    parser.add_argument('fasta', type=argparse.FileType("r"), help='Fasta file')
    parser.add_argument('--protein', action='store_true', help='The sequence is protein')
    parser.add_argument('--func', choices=funcs)
    args = parser.parse_args()
    main(**vars(args))
