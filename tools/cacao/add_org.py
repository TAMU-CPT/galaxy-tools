#!/usr/bin/env python
import json
import requests
import argparse
from Bio import SeqIO
from BCBio import GFF


def auth(creds, url):
    data = json.load(creds)['cacao']
    r = requests.post(url + 'api-token-auth/', data=data)
    return 'JWT ' + r.json()['token']


def get(token, url):
    q = requests.get(
        url,
        headers={'Authorization': token}
    ).json()
    return q


def post(token, url, data):
    q = requests.post(
        url,
        data=data,
        headers={'Authorization': token}
    ).json()
    return q


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('url', help='Backend URL. Include trailing slash.')
    parser.add_argument('creds', type=argparse.FileType("r"), help='json file with username/password')

    parser.add_argument('--taxon', type=str, default='NA')
    parser.add_argument('--ebi_id', type=str, default='NA')
    parser.add_argument('gff3', type=argparse.FileType('r'))
    parser.add_argument('fasta', type=argparse.FileType('r'))
    parser.add_argument('organism_file', type=argparse.FileType('r'))

    args = parser.parse_args()

    token = auth(args.creds, args.url)

    organism_name = json.load(args.organism_file)[0]
    if 'commonName' in organism_name:
        organism_name = organism_name['commonName']
    elif 'common_name' in organism_name:
        organism_name = organism_name['common_name']
    else:
        raise Exception("Bad input organism data")

    organism = post(token, args.url + 'organisms/', dict(
        common_name=organism_name,
        taxon=args.taxon,
        ebi_id=args.ebi_id
    ))
    if organism.get('common_name', [None])[0] == 'organism with this common name already exists.':
        organism = get(token, args.url + 'organisms/?common_name=' + organism_name)['results'][0]
        print("Organism[%s]: Pre-existing" % organism_name)
    else:
        print("Organism[%s]: Registered" % organism_name)

    refseqs = {}
    for record in SeqIO.parse(args.fasta, "fasta"):
        refseq = post(token, args.url + 'refseq/', dict(  # noqa
            name=record.id,
            length=len(record.seq),
            organism=organism['id'],
        ))
        if refseq.get('non_field_errors', [None])[0] == 'The fields name, organism must make a unique set.':
            refseq = get(token, args.url + 'refseq/?name=%s&%s' % (record.id, organism['id']))['results'][0]
        print("\tRefSeq[%s]: Pre-existing" % record.id)
    else:
        print("\tRefSeq[%s]: Registered" % record.id)

        refseqs[record.id] = refseq

    for rec in GFF.parse(args.gff3):
        # rs = RefSeq.objects.get(name=rec.id, organism=organism)
        rs = refseqs[rec.id]
        for feat in rec.features:
            if feat.type != 'gene':
                continue

            post(token, args.url + 'genes/', {
                "id": feat.id,
                "start": int(feat.location.start),
                "end": int(feat.location.end),
                "strand": int(feat.location.strand),
                "refseq": rs['id'],
                "db_object_id": feat.id,
                "db_object_symbol": feat.id,
                "db_object_name": "",
                "db_object_synonym": "",
                "db_object_type": "protein",
                "gene_product_id": "",
            })
            print("\tFeature[%s]: Updated" % feat.id)

if __name__ == '__main__':
    main()
