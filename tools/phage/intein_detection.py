#!/usr/bin/env
import sys
import argparse
import itertools
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


class InteinFinder(object):

    def __init__(self, genes, blast):
        feats = {}
        for rec in GFF.parse(genes):
            self.recid = rec.id
            self.reca = rec.annotations
            for feature in rec.features:
                feats[feature.id] = feature
        self.genes = feats
        self.blast = blast

    def _dist_between_genes(self, a, b):
        # http://stackoverflow.com/a/16843530/347368
        r1 = (a.location.start, a.location.end)
        r2 = (b.location.start, b.location.end)
        x, y = sorted((r1, r2))

        if x[0] <= x[1] < y[0] and all(y[0] <= y[1] for y in (r1, r2)):
            return y[0] - x[1]
        return 0

    def filter_possible_inteins(self, a_id, b_id, threshold=60):
        try:
            a = self.genes[a_id]
            b = self.genes[b_id]
        except KeyError:
            return False

        if a.location.strand != b.location.strand:
            return False

        if self._dist_between_genes(a, b) <= threshold:
            return True

        return False

    def intein_detection(self):
        data = {}
        for row in self.blast:
            rowdata = row.split('\t')
            gi_list = [rowdata[1]] + rowdata[12].split(';')
            for gi in gi_list:
                if gi in data:
                    if rowdata[0] not in data[gi]:
                        data[gi].append(rowdata[0])
                else:
                    data[gi] = [rowdata[0]]

        output = {}
        for k in data.keys():
            if len(data[k]) > 1:
                for x, y in itertools.combinations(data[k], 2):
                    if self.filter_possible_inteins(x, y):
                        ok = '%s..%s' % (x, y)
                        if ok in output:
                            output[ok]['hits'].append(k)
                        else:
                            output[ok] = {
                                'a': self.genes[x],
                                'b': self.genes[y],
                                'hits': [k]
                            }
        return output

    def output_to_gff3(self, id_output):
        seq = Seq("ACTG")
        rec = SeqRecord(seq, self.recid)

        for idx, possibility in enumerate(id_output):
            gene_a = id_output[possibility]['a']
            gene_b = id_output[possibility]['b']

            intein_start = min(gene_a.location.start, gene_b.location.start)
            intein_end = max(gene_a.location.end, gene_b.location.end)
            strand = gene_a.location.strand

            rec.features.append(
                SeqFeature(
                    FeatureLocation(intein_start, intein_end),
                    type="misc_feature",
                    strand=strand,
                    qualifiers={
                        "source": "InteinFinder",
                        "evidence": id_output[possibility]['hits'],
                        "ID": "intein_%s" % idx
                    }
                )
            )
        rec.annotations = self.reca

        return [rec]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Intein detection')
    parser.add_argument('genes', type=file, help='GFF3 Gene Calls')
    parser.add_argument('blast', type=file, help='Blast TSV Data')
    args = parser.parse_args()

    intf = InteinFinder(**vars(args))
    GFF.write(intf.output_to_gff3(intf.intein_detection()), sys.stdout)
