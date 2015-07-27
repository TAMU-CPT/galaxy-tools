import sys
import argparse
from BCBio import GFF


def rec_prefix_name(feature, prefix):
    # Stupid NCBI....
    if 'Derives_from' in feature.qualifiers:
        feature.qualifiers['Parent'] = [feature.qualifiers['Derives_from'][0].replace('.t', '.e')]
        del feature.qualifiers['Derives_from']

    # MORE stupid ncbi/genbank2gff3.pl
    if feature.type == 'exon':
        parent = feature.qualifiers['Parent'][0]
        feature.qualifiers['ID'] = [parent.replace('.t', '.e')]

    for attrib in ('Name', 'Parent', 'ID', 'Derives_from'):
        if attrib in feature.qualifiers:
            for idx, name in enumerate(feature.qualifiers[attrib]):
                feature.qualifiers[attrib][idx] = prefix + name

    if hasattr(feature, 'sub_features'):
        for subfeature in feature.sub_features:
            rec_prefix_name(subfeature, prefix)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Pre-process gff files')
    parser.add_argument('gff', type=file, help='gff3 file to process')
    parser.add_argument('--organism', type=str, help='organism')
    parser.add_argument('--source', type=str, help='adjust source value')
    parser.add_argument('--feature_prefix', type=str, help='adjust source value')
    args = parser.parse_args()

    for rec in GFF.parse(args.gff):
        original_name = rec.id
        if args.source is not None:
            rec.id = args.source
        else:
            rec.id = args.organism

        for feature in rec.features:
            if feature.type in ('region', 'chromosome'):
                feature.type = 'chromosome'
                feature.qualifiers['Name'] = [rec.id]
                feature.qualifiers['ID'] = [rec.id]

            if feature.type == 'polypeptide':
                feature.type = 'CDS'

            if args.feature_prefix is not None and not feature.type is 'chromosome':
                rec_prefix_name(feature, args.feature_prefix)


        GFF.write([rec], sys.stdout)
