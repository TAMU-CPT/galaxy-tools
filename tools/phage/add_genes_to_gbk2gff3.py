import sys
import argparse
from BCBio import GFF
from Bio.SeqFeature import SeqFeature, FeatureLocation
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def fix(gff):
    for record in GFF.parse(gff):
        feature_tag_groups = {}
        rejects = []
        for feature in record.features:
            fn = feature.qualifiers.get('Name', [None])[0]
            if fn is not None:
                if fn in feature_tag_groups:
                    feature_tag_groups[fn].append(feature)
                else:
                    feature_tag_groups[fn] = [feature]
            else:
                rejects.append(feature)

        for group in feature_tag_groups:
            if len(feature_tag_groups[group]) == 2:
                # GOOD! CDS + RBS
                starts = min([x.location.start for x in feature_tag_groups[group]])
                ends = max([x.location.end for x in feature_tag_groups[group]])
                cds = [x for x in feature_tag_groups[group] if x.type == 'CDS'][0]
                rbs = [x for x in feature_tag_groups[group] if x.type == 'RBS'][0]

                gene_feature = SeqFeature(
                    FeatureLocation(starts, ends, strand=cds.strand),
                    type="gene", id=cds.id, qualifiers={"ID": cds.id}
                )
                cds.qualifiers['ID'][0] += '.CDS'
                rbs.type = "Shine_Dalgarno_sequence"
                gene_feature.sub_features = [cds, rbs]

                rejects.append(gene_feature)
            else:
                log.warning("Found %s features in group %s",
                            len(feature_tag_groups[group]), group)
                rejects.extend(feature_tag_groups[group])

        record.features = rejects
        yield [record]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix converted files')
    parser.add_argument('gff', type=argparse.FileType("r"), help='bp_gbk2gff3 output')
    args = parser.parse_args()

    for record in fix(**vars(args)):
        record.annotations = {}
        GFF.write(record, sys.stdout)
