#!/usr/bin/env python
import BIO_FIX_TOPO  # NOQA
import sys
import random
import argparse
import logging
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()
random.seed(42)


def get_id(feature=None, parent_prefix=None):
    result = ""
    if parent_prefix is not None:
        result += parent_prefix + "|"

    for x in ("locus_tag", "protein_id", "gene"):
        if x in feature.qualifiers:
            result += feature.qualifiers[x][0]
            break
    else:
        return feature.id
    return result


def same_end(query, target):
    results = []
    for x in target:
        if x.location.strand == query.location.strand:
            if x.location.strand > 0 and x.location.end == query.location.end:
                results.append(x)
            elif x.location.strand < 0 and x.location.start == query.location.start:
                results.append(x)
    return results


def sortOrder(feature):
    # if 'locus_tag' in feature.qualifiers:
    # so = feature.qualifiers['locus_tag'][0] + '__' + "%010d" % feature.location.start
    # else:
    # so = "%010d" % feature.location.start

    return feature.location.start


def nearbyRbss(cds, rbss):
    for r in rbss:
        if cds.strand > 0:
            return [
                r
                for r in rbss
                if abs(r.location.end - cds.location.start) < 20
                and r.location.strand == cds.location.strand
            ]
        else:
            return [
                r
                for r in rbss
                if abs(r.location.start - cds.location.end) < 20
                and r.location.strand == cds.location.strand
            ]


def unionLoc(a, b):
    return FeatureLocation(
        start=min(a.start, b.start), end=max(a.end, b.end), strand=a.strand
    )


def fix_locus(g, cds):
    gene_locus = g.qualifiers.get("locus_tag", [None])[0]
    cds_locus = cds.qualifiers.get("locus_tag", [None])[0]
    if gene_locus and cds_locus:
        if gene_locus == cds_locus:
            pass
        else:
            cds.qualifiers["cpt_gmc"] = ["Different locus tag from associated gene."]
            cds.qualifiers["cpt_gmc_locus"] = cds.qualifiers["locus_tag"]
            cds.qualifiers["locus_tag"] = [gene_locus]
    elif gene_locus and not cds_locus:
        cds.qualifiers["cpt_gmc"] = ["Missing Locus Tag"]
        cds.qualifiers["locus_tag"] = [gene_locus]
    elif not gene_locus and cds_locus:
        g.qualifiers["cpt_gmc"] = ["Missing Locus Tag"]
        g.qualifiers["locus_tag"] = [cds_locus]
    else:
        locus_tag = "cpt_orf_%s" % random.randint(0, 100000)
        g.qualifiers["cpt_gmc"] = ["BOTH Missing Locus Tag"]
        g.qualifiers["locus_tag"] = [locus_tag]
        cds.qualifiers["cpt_gmc"] = ["BOTH Missing Locus Tag"]
        cds.qualifiers["locus_tag"] = [locus_tag]


def correct_model(genbank_file):
    for record in SeqIO.parse(genbank_file, "genbank"):
        cds = [f for f in record.features if f.type == "CDS"]
        trna = [f for f in record.features if f.type == "tRNA"]
        genes = [f for f in record.features if f.type == "gene"]
        rbss = [f for f in record.features if f.type == "RBS"]

        # No genes at all
        if len(genes) == 0:
            for c in cds:
                quals = {
                    "cpt_source": ["CPT_GENE_MODEL_CORRECTION"],
                    "gene": c.qualifiers.get("gene", []),
                    "product": c.qualifiers.get("product", []),
                    "locus_tag": c.qualifiers.get("locus_tag", [get_id(c)]),
                }
                if "locus_tag" not in c.qualifiers:
                    c.qualifiers["locus_tag"] = [get_id(c)]

                gene_location = None
                if len(rbss) > 0:
                    r = nearbyRbss(c, rbss)
                    if len(r) == 1:
                        gene_location = unionLoc(c.location, r[0].location)
                        r[0].qualifiers["locus_tag"] = c.qualifiers["locus_tag"]
                    else:
                        gene_location = FeatureLocation(
                            c.location.start, c.location.end, c.location.strand
                        )
                else:
                    gene_location = FeatureLocation(
                        c.location.start, c.location.end, c.location.strand
                    )

                record.features.append(
                    SeqFeature(location=gene_location, type="gene", qualifiers=quals)
                )

            # print genes
            record.features += genes
        elif len(genes) != len(cds):
            for g in genes:
                associated_cds = same_end(g, cds)
                if len(associated_cds) == 1:
                    # Ensure matching locus tags
                    fix_locus(g, associated_cds[0])
                elif len(associated_cds) == 0:
                    associated_trna = same_end(g, trna)
                    if len(associated_trna) == 1:
                        pass
                    else:
                        log.warn(
                            "Could not find a child feature for gene %s", get_id(g)
                        )
                else:
                    log.warn("Could not find a child feature for gene %s", get_id(g))

            log.info(
                "Different number of CDSs and genes. There may be bugs in this process (genes=%s != cds=%s)",
                len(genes),
                len(cds),
            )
        elif len(genes) == len(cds):
            for g in genes:
                associated_cds = same_end(g, cds)
                if len(associated_cds) == 1:
                    # Ensure matching locus tags
                    fix_locus(g, associated_cds[0])
                elif len(associated_cds) == 0:
                    associated_trna = same_end(g, trna)
                    if len(associated_trna) == 1:
                        pass
                    else:
                        log.warn(
                            "Could not find a child feature for gene %s", get_id(g)
                        )
                else:
                    log.warn("Could not find a child feature for gene %s", get_id(g))

        record.features = sorted(record.features, key=lambda x: sortOrder(x))
        yield [record]


if __name__ == "__main__":
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(
        description="Export a subset of features from a Genbank file", epilog=""
    )
    parser.add_argument(
        "genbank_file", type=argparse.FileType("r"), help="Genbank file"
    )

    args = vars(parser.parse_args())
    for seq in correct_model(**args):
        SeqIO.write(seq, sys.stdout, "genbank")
