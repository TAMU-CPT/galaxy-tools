#!/usr/bin/env python
import sys
import argparse
import logging
import uuid
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gff3 import feature_lambda, feature_test_type, get_id
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def main(fasta, gff3, feature_filter=None, nodesc=False):

    if feature_filter == 'nice_cds':
        from gff2gb import gff3_to_genbank
        
        for rec in gff3_to_genbank(gff3, fasta):
            seenList = {}
            if rec.seq[0] == '?':
              print("No Fasta ID matches GFF")
              exit(1) 
            for feat in sorted(rec.features, key=lambda x: x.location.start):
                if feat.type != 'CDS':
                    continue
                  
                ind = 0              
                if str(feat.qualifiers.get('locus_tag', get_id(feat)).replace(' ', '-')) in seenList.keys():
                  seenList[str(feat.qualifiers.get('locus_tag', get_id(feat)).replace(' ', '-'))] += 1
                  ind = seenList[str(feat.qualifiers.get('locus_tag', get_id(feat)).replace(' ', '-'))]
                else:
                  seenList[str(feat.qualifiers.get('locus_tag', get_id(feat)).replace(' ', '-'))] = 1
                append = ""
                if ind != 0:
                  append = "_" + str(ind)

                if nodesc:
                    description = ''
                else:
                    feat.qualifiers['ID'] = [feat._ID]
                    product = feat.qualifiers.get('product', '')
                    description = '{1} [Location={0.location};ID={0.qualifiers[ID][0]}]'.format(feat, product)
                #print(feat.qualifiers.get('locus_tag', get_id(feat)).replace(' ', '-'))
                yield [
                    SeqRecord(
                        feat.extract(rec).seq,
                        id=str(feat.qualifiers.get('locus_tag', get_id(feat)).replace(' ', '-')) + append,
                        description=description
                    )
                ]
                

    elif feature_filter == 'unique_cds':
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        seen_ids = {}
         
        for rec in GFF.parse(gff3, base_dict=seq_dict):
            noMatch = True
            for x in seq_dict:
              if x == rec.id:
                noMatch = False
            if noMatch:
              print("No Fasta ID matches GFF")
              exit(1)
            newfeats = []
            for feat in sorted(feature_lambda(
                rec.features,
                feature_test_type,
                {'type': 'CDS'},
                subfeatures=False
            ), key=lambda f: f.location.start):
                nid = rec.id + '____' + feat.id
                if nid in seen_ids:
                    nid = nid + '__' + uuid.uuid4().hex
                feat.qualifiers['ID'] = nid
                newfeats.append(feat)
                seen_ids[nid] = True

                if nodesc:
                    description = ''
                else:
                    important_data = {
                        'Location': feat.location,
                    }
                    if 'Name' in feat.qualifiers:
                        important_data['Name'] = feat.qualifiers.get('Name', [''])[0]

                    description = '[{}]'.format(
                        ';'.join([
                            '{key}={value}'.format(key=k, value=v) for (k, v) in important_data.items()
                        ])
                    )

                yield [
                    SeqRecord(
                        feat.extract(rec).seq,
                        id=nid.replace(' ', '-'),
                        description=description
                    )
                ]
            rec.features = newfeats
            rec.annotations = {}
            GFF.write([rec], sys.stderr)
    else:
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        for rec in GFF.parse(gff3, base_dict=seq_dict):
            noMatch = True
            for x in seq_dict:
              if x == rec.id:
                noMatch = False
            if noMatch:
              print("No Fasta ID matches GFF")
              exit(1)
            for feat in sorted(feature_lambda(
                rec.features,
                feature_test_type,
                {'type': feature_filter},
                subfeatures=False
            ), key=lambda f: f.location.start):
                id = feat.id
                if len(id) == 0:
                    id = get_id(feat)

                if nodesc:
                    description = ''
                else:
                    important_data = {
                        'Location': feat.location,
                    }
                    if 'Name' in feat.qualifiers:
                        important_data['Name'] = feat.qualifiers.get('Name', [''])[0]

                    description = '[{}]'.format(
                        ';'.join([
                            '{key}={value}'.format(key=k, value=v) for (k, v) in important_data.items()
                        ])
                    )

                yield [
                    SeqRecord(
                        feat.extract(rec).seq,
                        id=id.replace(' ', '-'),
                        description=description
                    )
                ]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export corresponding sequence in genome from GFF3', epilog="")
    parser.add_argument('fasta', type=argparse.FileType("r"), help='Fasta Genome')
    parser.add_argument('gff3', help='GFF3 File')
    parser.add_argument('--feature_filter', default=None, help='Filter for specific feature types')
    parser.add_argument('--nodesc', action='store_true', help='Strip description field off')
    args = parser.parse_args()

    for seq in main(**vars(args)):
        SeqIO.write(seq, sys.stdout, 'fasta')
