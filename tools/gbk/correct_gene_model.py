#!/usr/bin/env python
import sys
import argparse
import logging
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def get_id(feature=None, parent_prefix=None):
    result = ""
    if parent_prefix is not None:
        result += parent_prefix + '|'

    for x in ('locus_tag', 'protein_id', 'gene'):
        if x in feature.qualifiers:
            result += feature.qualifiers[x][0]
            break
    else:
        return feature.id
    return result


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
        trna = [f for f in record.features if f.type == 'tRNA']
        genes = [f for f in record.features if f.type == 'gene']
        rbss = [f for f in record.features if f.type == 'RBS']

        # No genes at all
        if len(genes) == 0:
            for c in cds:
                quals = {
                    'cpt_source': ['CPT_GENE_MODEL_CORRECTION'],
                    'gene': c.qualifiers.get('gene', []),
                    'product': c.qualifiers.get('product', []),
                    'locus_tag': c.qualifiers.get('locus_tag', [get_id(c)]),
                    'ID': get_id(c),
                }
                if 'gene' not in c.qualifiers:
                    c.qualifiers['gene'] = quals['ID']
                if 'locus_tag' not in c.qualifiers:
                    c.qualifiers['locus_tag'] = get_id(c)

                if len(rbss) > 0:
                    extra_feats = []
                    r = nearbyRbss(c, rbss)
                    if len(r) == 1:
                        extra_feats.append(SeqFeature(
                            unionLoc(c.location, r[0].location),
                            type="gene",
                            qualifiers=quals
                        ))
                    else:
                        extra_feats.append(SeqFeature(
                            c.location,
                            type="gene",
                            qualifiers=quals
                        ))
                    record.features.extend(extra_feats)
                else:
                    record.features.append(
                        SeqFeature(
                            c.location,
                            type="gene",
                            strand=c.location.strand,
                            qualifiers=quals
                        )
                    )

            # print genes
            record.features += genes
        elif len(genes) != len(cds):
            for g in genes:
                associated_cds = [x for x in cds if x.location.end == g.location.end]
                if len(associated_cds) == 1:
                    pass
                elif len(associated_cds) == 0:
                    associated_trna = [x for x in trna if x.location.end == g.location.end]
                    if len(associated_trna) == 1:
                        pass
                    else:
                        log.warn("Could not find a child feature for gene %s", get_id(g))
                else:
                    log.warn("Could not find a child feature for gene %s", get_id(g))

            log.info("Different number of CDSs and genes. There may be bugs in this process (genes=%s != cds=%s)", len(genes), len(cds))

        record.features = sorted(record.features, key=lambda x: int(x.location.start) - (1 if x.type == 'gene' else 0))
        yield [record]


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Export a subset of features from a Genbank file', epilog="")
    parser.add_argument('genbank_file', type=file, help='Genbank file')

    args = vars(parser.parse_args())
    for seq in correct_model(**args):
        SeqIO.write(seq, sys.stdout, 'genbank')
