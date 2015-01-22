#!/usr/bin/env python
import argparse
import StringIO
import logging
logging.basicConfig(level=logging.INFO)
from BCBio import GFF
from Bio import SeqIO
import hashlib

def uniqify(some_list):
    return sorted(list(set(some_list)))

def csl(comma_separated_list):
    """
    This is overkill
    """
    if len(comma_separated_list) == 1:
        return comma_separated_list[0].split('","')
    elif len(comma_separated_list) > 1:
        # recurse
        results = []
        for item in comma_separated_list:
            results.extend(comma_separated_list(item))
        return results
    else:
        return []

def extract_features(gff_file=None, **kwd):
    if gff_file is None:
        raise ValueError("Must specify gff file")

    feature_qualifiers = {}

    for rec in GFF.parse(gff_file):
        quals = {
            'dbxref': [],
            'note': [],
        }
        md5 = None
        for feature in rec.features:
            # Ehhhh...assume that NCBI wants IPR numbers and not other databases
            #if 'Name' in feature.qualifiers:
                #quals['dbxref'].extend(feature.qualifiers['Name'])

            if 'Ontology_term' in feature.qualifiers:
                quals['dbxref'].extend(csl(feature.qualifiers['Ontology_term']))

            if 'Dbxref' in feature.qualifiers:
                quals['dbxref'].extend(csl(feature.qualifiers['Dbxref']))

            quals['id'] = feature.id

            # Might be useful as a note?
            #if 'signature_desc' in feature.qualifiers:
                #quals['note'].extend(feature.qualifiers['signature_desc'])

            if 'md5' in feature.qualifiers:
                md5 = feature.qualifiers['md5'][0]

        if md5 is not None:
            quals['note'] = uniqify(quals['note'])
            quals['dbxref'] = uniqify(quals['dbxref'])
            feature_qualifiers[md5] = quals
    return feature_qualifiers


def hash(sequence):
    return str(hashlib.md5(str(sequence)).hexdigest())

def merge_features(ipr_features=None, genbank_file=None, **kwd):
    output = StringIO.StringIO()
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            feature_hash = hash(feature.extract(record.seq).translate(table=11).strip('*'))
            #print feature.extract(record.seq).translate(table=11).strip('*'), feature_hash
            if feature_hash in ipr_features:
                if 'dbxref' not in feature.qualifiers:
                    feature.qualifiers['dbxref'] = []
                if 'noteg' not in feature.qualifiers:
                    feature.qualifiers['note'] = []

                feature.qualifiers['note'] += ipr_features[feature_hash]['note']
                feature.qualifiers['dbxref'] += ipr_features[feature_hash]['dbxref']
        SeqIO.write(record, output, "genbank")
    return output.getvalue()


__doc__ = """
Merge InterProScan results into a GenBank File
==============================================

"""

if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Merge InterProScan data into a Genbank file')
    parser.add_argument('genbank_file', type=file, help='Genbank file')
    parser.add_argument('gff_file', type=file, help='InterPro GFF3 Input')

    args = vars(parser.parse_args())

    features = extract_features(**args)
    result = merge_features(ipr_features=features, **args)
    print result
