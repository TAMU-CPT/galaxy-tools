#!/usr/bin/env python
"""Blast based reopening class"""
import sys
import subprocess
import numpy as np
from Bio import SeqIO


class BlastBasedReopening(object):
    """Re-open genomes based on blast"""

    def __init__(self, genome="reopen-data/mis.fa"):
        self.genome = genome
        genome_seq = SeqIO.read(self.genome, "fasta")

        blast_results = self.getBlast(self.genome)
        # If more than half the hits are 'backwards'

        if self.shouldReverse(blast_results) > 0.5:
            new_genome_file = self.genome.replace(".fa", ".revcom.fa")
            SeqIO.write(
                [genome_seq.reverse_complement(id=True, description=True)],
                new_genome_file,
                "fasta",
            )
            blast_results = self.getBlast(new_genome_file)

        ref_genome_hits = [
            self.avg(hit["sstart"], hit["send"]) for hit in blast_results
        ]
        our_genome_hits = [
            self.avg(hit["qstart"], hit["qend"]) for hit in blast_results
        ]
        # Now we sort relative to ref genome
        hits = zip(our_genome_hits, ref_genome_hits)
        hits = sorted(hits, key=lambda x: x[1])
        # Then we re-construct *_genome_hits from the sorted/zipped version
        ref_genome_hits = [x[1] for x in hits]
        our_genome_hits = [x[0] for x in hits]
        # These are now two independent lists, sorted relative to ref genome. We
        # can shift our_genome_hits without affecting ref genome ordering.

        (score, (index, m, c)) = self.scoreResults(ref_genome_hits, our_genome_hits)

        # Now on to plotting the real data.
        totalData = [
            (hit["qstart"], hit["qend"], hit["sstart"], hit["send"])
            for hit in blast_results
        ]
        totalData = sorted(totalData, key=lambda x: self.avg(x[2], x[3]))
        print "Reopen at", totalData[0][0]

    def getBlast(self, path):
        blast_results = (
            subprocess.check_output(
                [
                    "blastn",
                    "-db",
                    "test-data/canonical_nucl",
                    "-query",
                    path,
                    "-outfmt",
                    "6 sseqid evalue pident qstart qend sstart send",
                ]
            )
            .strip()
            .split("\n")
        )
        blast_results = [x.split("\t") for x in blast_results]
        blast_results = [
            {
                "sseqid": sseqid,
                "evalue": float(evalue),
                "pident": float(pident),
                "qstart": int(qstart),
                "qend": int(qend),
                "sstart": int(sstart),
                "send": int(send),
            }
            for (sseqid, evalue, pident, qstart, qend, sstart, send) in blast_results
        ]
        if len([x["sseqid"] for x in blast_results]) > 1:
            sys.stderr.write("DO not know how to handle. This WILL MISBEHAVE\n")
        return blast_results

    def shouldReverse(self, blast_results):
        should_reverse = 0
        total = len(blast_results)

        for hit in blast_results:
            qdir = 1 if hit["qstart"] > hit["qend"] else -1
            sdir = 1 if hit["sstart"] > hit["send"] else -1
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
