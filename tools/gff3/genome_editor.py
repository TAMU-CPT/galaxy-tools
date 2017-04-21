#!/usr/bin/env python
import logging
import copy
import argparse
import tsv
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def mutate(gff3, fasta, changes, new_id):
    # Change Language
    # - we can only accept ONE genome as an input. (TODO: support multiple?)
    # - we can only build ONE genome as an output. (TODO: support multiple?)
    # - must allow selection of various regions
    # '1,1000,+;40,100,-'
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    # Pull first and onl record
    rec = list(GFF.parse(gff3, base_dict=seq_dict))[0]
    # Create a "clean" record
    new_record = copy.deepcopy(rec)
    new_record.id = new_id
    new_record.seq = Seq('')
    new_record.features = []
    new_record.annotations = {}
    # Process changes.
    chain = []
    for change in changes:
        (start, end, strand) = change.split(',')
        start = int(start) - 1
        end = int(end)
        if strand == '+':
            tmp_req = rec[start:end]
        else:
            tmp_req = rec[start:end].reverse_complement(id=True, name=True, description=True, features=True, annotations=True, letter_annotations=True, dbxrefs=True)

        chain.append([
            rec.id,
            start,
            end,
            strand,
            new_record.id,
            len(new_record),
            len(new_record) + (end - start),
            '+'
        ])

        new_record.seq += tmp_req.seq
        # NB: THIS MUST USE BIOPYTHON 1.67. 1.68 Removes access to
        # subfeatures, which means you will only get top-level features.
        new_record.features += tmp_req.features
    yield new_record, chain

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gff3', type=argparse.FileType("r"), help='Annotations')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='Sequence')
    parser.add_argument('new_id', help='Append to ID', default='_v2')
    parser.add_argument('--out_fasta', type=argparse.FileType("w"), help='Output fasta', default='out.fa')
    parser.add_argument('--out_gff3', type=argparse.FileType("w"), help='Output gff3', default='out.gff3')
    parser.add_argument('--out_simpleChain', type=argparse.FileType("w"), help='Output simple chain (i.e. not a real UCSC chain file)', default='out.chain')
    parser.add_argument('--changes', nargs='+')
    args = parser.parse_args()

    for rec, chain in mutate(args.gff3, args.fasta, args.changes, args.new_id):
        # TODO: Check that this appends and doesn't overwirte
        GFF.write([rec], args.out_gff3)
        SeqIO.write([rec], args.out_fasta, 'fasta')
        tsv.dump(chain, args.out_simpleChain)
