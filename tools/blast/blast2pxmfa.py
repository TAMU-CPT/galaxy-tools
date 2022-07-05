#!/usr/bin/env python
import os
import argparse
from CPT_GFFParser import gffParse, gffWrite
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gff3 import feature_lambda, feature_test_true
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import tempfile
import logging

log = logging.getLogger()
log.setLevel(logging.DEBUG)


def parse_gff3(annotations, genome):
    annotations.seek(0)
    genome.seek(0)

    data = {}
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    for record in gffParse(annotations, base_dict=seq_dict):
        for feature in feature_lambda(
            record.features, feature_test_true, {}, subfeatures=False
        ):
            data[feature.id] = {
                "rec": record.id,
                "loc": feature.location,
                "seq": feature.extract(record).seq.translate(table=11, cds=False)[0:-1],
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
        line = row.split("\t")
        yield {"from": line[0], "to": line[1], "pident": line[2], "evalue": line[10]}


def cluster_relationships(data):
    generated_clusters_ids = []

    for relationship in data:
        hit = False
        relationship_set = set((relationship["from"], relationship["to"]))
        for idx, cluster in enumerate(generated_clusters_ids):
            if relationship["from"] in cluster or relationship["to"] in cluster:
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

    # Clustal mangles IDs. Fricking Clustal. Specifically it truncates them
    # meaning we have to RE-ID every single sequence that goes into it. Lovely.
    id_map = {}
    for idx, seq in enumerate(SequenceList):
        id_map[str(idx)] = seq.id
        seq.id = str(idx)

    t_fa = tempfile.NamedTemporaryFile(prefix="blastxmfa.", delete=False)
    t_aln = tempfile.NamedTemporaryFile(prefix="blastxmfa.", delete=False)
    # Write fasta to file
    SeqIO.write(SequenceList, str(t_fa.name), "fasta")
    # Flush & close
    t_fa.flush()
    t_fa.close()
    t_aln.close()

    # Build clustal CLI
    d = ClustalwCommandline(infile=t_fa.name, outfile=t_aln.name)
    logging.debug(d)
    # Call clustal
    d()
    # Get our alignment back
    try:
        aln = AlignIO.read(t_aln.name, "clustal")
        # Now we replace the IDs with the correct, full length ones
        for a in aln:
            a.id = id_map[a.id]
        # Cleanup
        os.unlink(t_fa.name)
        os.unlink(t_aln.name)
        return aln
    except Exception as e:
        logging.error("%s, %s", e, t_fa.name)
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


def smaller_than_3n(x, it):
    THRESH = 3 * x
    for item in it:
        if len(item) <= THRESH:
            yield item
        else:
            log.warn(
                "Cluster with %s (>%s) members, seems excessive", len(item), THRESH
            )


class XmfaWriter(object):
    HEADER_TPL = "> {seqnum}:{start}-{end} {strand} {file} # {realname}\n"

    def __init__(self, handle):
        self.output = handle
        self.output.write("#FormatVersion Mauve1\n")

    def write(self, seqnum, start, end, strand, file, realname, sequence):
        self.output.write(
            self.HEADER_TPL.format(
                seqnum=seqnum,
                start=start,
                end=end,
                strand=strand,
                file=file,
                realname=realname,
            )
        )

        for line in split_by_n(sequence, 80):
            self.output.write(line + "\n")

    def end(self):
        self.output.write("=\n")


def blast2pxmfa(blast, fasta, gff3, output, genomic=False):
    logging.info("Parsing sequence")
    locations = parse_gff3(gff3, fasta)
    logging.info("Parsed locations, clustering")
    recids = get_seqids(fasta)
    logging.info("Found %s records in fasta", len(recids))

    xmw = XmfaWriter(output)

    # First let's generate some blastclust style clusters.
    clusters = list(
        smaller_than_3n(
            len(recids),
            larger_than_one(cluster_relationships(gen_relationships(blast))),
        )
    )
    logging.debug("%s clusters generated", len(clusters))

    def sortIndexForFeatId(element):
        # Will not work correctly if genome length is >1mb
        general = recids.index(locations[element]["rec"]) * 1000000
        specific = locations[element]["loc"].start
        return general + specific

    for idx, cluster in enumerate(clusters):
        logging.debug("Cluster %s/%s, size=%s", idx + 1, len(clusters), len(cluster))
        # We're considering 1 LCB :: 1 cluster
        seqs = []
        for element in cluster:
            if element not in locations:
                logging.warning("Could not find this feature %s", element)
                continue

            sr = SeqRecord(
                locations[element]["seq"],
                id=element,
                description="[{0.start}:{0.end}:{0.strand}]".format(
                    locations[element]["loc"]
                ),
            )
            seqs.append(sr)

        aligned_seqs = align_sequences(seqs)
        if aligned_seqs is None:
            logging.error("Error aligning cluster [%s]", "".join(cluster))
            continue

        sortedCluster = sorted(cluster, key=lambda x: sortIndexForFeatId(x))
        sortedAligned = sorted(aligned_seqs, key=lambda x: sortIndexForFeatId(x.id))
        # print(sortedCluster, [x.id for x in sortedAligned])

        # Pre-check the asserts.
        goodToGo = all(
            [
                element == aligned_seq.id
                for (element, aligned_seq) in zip(sortedCluster, sortedAligned)
            ]
        )
        if not goodToGo:
            logging.info(
                "Skipping one grouping: %s != %s",
                ",".join(sortedCluster),
                ",".join([x.id for x in sortedAligned]),
            )
            # Skip this look
            continue
        # print(aligned_seqs)
        for element, aligned_seq in zip(sortedCluster, sortedAligned):
            if element not in locations:
                logging.warning("Could not find this feature %s", element)
                continue

            eloc = locations[element]
            xmw.write(
                seqnum=recids.index(eloc["rec"]) + 1,
                start=eloc["loc"].start,
                end=eloc["loc"].end,
                strand="+" if eloc["loc"].strand == 1 else "-",
                file=eloc["rec"] + ".fa",
                realname=element,
                sequence=str(aligned_seq.seq),
            )
        xmw.end()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert Blast TSV to protein XMFA", epilog=""
    )
    parser.add_argument("blast", type=argparse.FileType("r"), help="Blast TSV Output")
    parser.add_argument("fasta", type=argparse.FileType("r"), help="Blast Input Fasta")
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 Gene Calls")
    parser.add_argument(
        "--genomic",
        action="store_true",
        help="Further reduce protein results into genomic-level results.",
    )
    parser.add_argument(
        "output", type=argparse.FileType("w"), help="Output file or - for stdout"
    )
    args = parser.parse_args()

    blast2pxmfa(**vars(args))
