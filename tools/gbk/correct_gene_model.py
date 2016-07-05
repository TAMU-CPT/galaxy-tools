#!/usr/bin/env python
import sys
import argparse
import logging
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def nearbyRbss(cds, rbss):
    for r in rbss:
        if cds.strand > 0:
            return [r for r in rbss if
                    abs(r.location.end - cds.location.start) < 20 and r.location.strand == cds.location.strand]
        else:
            return [r for r in rbss if
                    abs(r.location.start - cds.location.end) < 20 and r.location.strand == cds.location.strand]


def unionLoc(a, b):
    return FeatureLocation(
        start=min(a.start, b.start),
        end=max(a.end, b.end),
        strand=a.strand
    )


def correct_model(genbank_file):
    for record in SeqIO.parse(genbank_file, "genbank"):
        cds = [f for f in record.features if f.type == 'CDS']
        genes = [f for f in record.features if f.type == 'gene']
        rbss = [f for f in record.features if f.type == 'RBS']

        # No genes at all
        if len(genes) == 0:
            if len(rbss) > 0:
                extra_feats = []
                for c in cds:
                    r = nearbyRbss(c, rbss)
                    print r
                    if len(r) == 1:
                        extra_feats.append(SeqFeature(
                            unionLoc(c.location, r[0].location),
                            type="gene",
                            qualifiers={
                                'gene': f.qualifiers.get('gene', []),
                                'product': f.qualifiers.get('product', []),
                                'protein_id': f.qualifiers.get('protein_id', []),
                                'locus_tag': f.qualifiers.get('locus_tag', []),
                            }
                        ))
                    else:
                        extra_feats.append(SeqFeature(
                            c.location,
                            type="gene",
                            qualifiers={
                                'gene': f.qualifiers.get('gene', []),
                                'product': f.qualifiers.get('product', []),
                                'protein_id': f.qualifiers.get('protein_id', []),
                                'locus_tag': f.qualifiers.get('locus_tag', []),
                            }
                        ))
                record.features.extend(extra_feats)
            else:
                genes = [
                    SeqFeature(
                        f.location,
                        type="gene",
                        strand=f.location.strand,
                        qualifiers={
                            'gene': f.qualifiers.get('gene', []),
                            'product': f.qualifiers.get('product', []),
                            'protein_id': f.qualifiers.get('protein_id', []),
                            'locus_tag': f.qualifiers.get('locus_tag', []),
                        }
                    )
                    for f in cds
                ]
            # print genes
            record.features += genes
        elif len(genes) != len(cds):
            log.error("Different number of CDSs and genes. I don't know how to handle this case (yet)")

        record.features = sorted(record.features, key=lambda x: x.location.start)
        yield [record]


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Export a subset of features from a Genbank file', epilog="")
    parser.add_argument('genbank_file', type=file, help='Genbank file')

    args = vars(parser.parse_args())
    for seq in correct_model(**args):
        SeqIO.write(seq, sys.stdout, 'genbank')
