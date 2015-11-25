#!/usr/bin/env python
import os
import sys
import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gff3 import feature_lambda, feature_test_true
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import tempfile
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name='blast2xmfap')


HEADER_TPL = '> {seqnum}:{start}-{end} {strand} {file} # {realname}'


def parse_gff3(annotations, genome):
    annotations.seek(0)
    genome.seek(0)

    data = {}
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    for record in GFF.parse(annotations, base_dict=seq_dict):
        for feature in feature_lambda(
            record.features,
            feature_test_true,
            {},
            subfeatures=False
        ):
            data[feature.id] = {
                'rec': record.id,
                'loc': feature.location,
                'seq': feature.extract(record).seq.translate(table=11)[0:-1],
            }
    return data


def get_seqids(genome):
    genome.seek(0)
    recids = []
    for record in SeqIO.parse(genome, "fasta"):
        recids.append(record.id)
    return recids


def gen_relationships(blast):
    for row in blast:
        line = row.split('\t')
        yield {
            'from': line[0],
            'to': line[1],
            'pident': line[2],
            'evalue': line[10],
        }


def cluster_relationships(data):
    generated_clusters_ids = []

    for relationship in data:
        hit = False
        relationship_set = set((relationship['from'], relationship['to']))
        for idx, cluster in enumerate(generated_clusters_ids):
            if relationship['from'] in cluster or relationship['to'] in cluster:
                generated_clusters_ids[idx] = cluster.__or__(relationship_set)
                hit = True
                break

        if not hit:
            generated_clusters_ids.append(relationship_set)
    return generated_clusters_ids


def align_sequences(SequenceList):
    # No sense in aligning one sequence.
    if len(SequenceList) < 2:
        return None

    t_fa = tempfile.NamedTemporaryFile(prefix='blastxmfa.', delete=False)
    t_aln = tempfile.NamedTemporaryFile(prefix='blastxmfa.', delete=False)
    # Write fasta to file
    SeqIO.write(SequenceList, t_fa, 'fasta')
    # Flush & close
    t_fa.flush()
    t_fa.close()
    t_aln.close()

    # Build clustal CLI
    d = ClustalwCommandline(infile=t_fa.name, outfile=t_aln.name)
    # Call clustal
    d()
    # Get our alignment back
    try:
        aln = AlignIO.read(t_aln.name, 'clustal')
        # Cleanup
        os.unlink(t_fa.name)
        os.unlink(t_aln.name)
        return aln
    except Exception, e:
        log.error("%s, %s", e, t_fa.name)
        os.unlink(t_aln.name)
        return None


def split_by_n(seq, n):
    """A generator to divide a sequence into chunks of n units."""
    # http://stackoverflow.com/questions/9475241/split-python-string-every-nth-character
    while seq:
        yield seq[:n]
        seq = seq[n:]


def larger_than_one(it):
    for item in it:
        if len(item) > 1:
            yield item


def blast2xmfap(blast, fasta, gff3):
    locations = parse_gff3(gff3, fasta)
    log.debug("Parsed locations, clustering")
    recids = sorted(get_seqids(fasta))
    log.debug("Found %s records in fasta", len(recids))

    # First let's generate some blastclust style clusters.
    clusters = list(larger_than_one(cluster_relationships(gen_relationships(blast))))
    log.debug("%s clusters generated", len(clusters))

    print "#FormatVersion Mauve1"
    for cluster in clusters:
        # We're considering 1 LCB :: 1 cluster
        seqs = []
        for element in cluster:
            if element not in locations:
                log.warning("Could not find this feature %s", element)
                continue

            sr = SeqRecord(
                locations[element]['seq'],
                id=element,
                description='[{0.start}:{0.end}:{0.strand}]'.format(locations[element]['loc'])
            )
            seqs.append(sr)

        aligned_seqs = align_sequences(seqs)
        if aligned_seqs is None:
            log.error("Error aligning cluster [%s]", ''.join(cluster))
            continue

        for element, aligned_seq in zip(sorted(cluster), sorted(aligned_seqs, key=lambda x: x.id)):
            # Must be the same or we have big problems
            assert element == aligned_seq.id

            if element not in locations:
                log.warning("Could not find this feature %s", element)
                continue


            eloc = locations[element]
            print HEADER_TPL.format(
                seqnum=recids.index(eloc['rec']) + 1,
                start=eloc['loc'].start,
                end=eloc['loc'].end,
                strand='+' if eloc['loc'].strand == 1 else '-',
                file=eloc['rec'] + '.fa',
                realname=element,
            )

            for line in split_by_n(str(aligned_seq.seq), 80):
                print line

        print '='



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Blast 25 column TSV to XMFAp', epilog='')
    parser.add_argument('blast', type=file, help='Blast TSV Output')
    parser.add_argument('fasta', type=file, help='Blast Input Fasta')
    parser.add_argument('gff3', type=file, help='GFF3 Gene Calls')
    args = parser.parse_args()

    blast2xmfap(**vars(args))
