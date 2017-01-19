import json
import requests
import argparse
from Bio import SeqIO
from BCBio import GFF


def auth(creds, url):
    data = json.load(creds)['cacao']
    r = requests.post(url + 'api-token-auth/', data=data)
    return 'JWT ' + r.json()['token']


def post(token, url, data):
    q = requests.post(
        url,
        data=data,
        headers={'Authorization': token}
    ).json()
    print(q)
    return q


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('url', help='Backend URL. Include trailing slash.')
    parser.add_argument('creds', type=argparse.FileType("r"), help='json file with username/password')

    parser.add_argument('--organism_name', type=str)
    parser.add_argument('--taxon', type=str, default='NA')
    parser.add_argument('--ebi_id', type=str, default='NA')
    parser.add_argument('gff3', type=argparse.FileType('r'))
    parser.add_argument('fasta', type=argparse.FileType('r'))

    args = parser.parse_args()

    token = auth(args.creds, args.url)

    organism = post(token, args.url + 'organisms/', dict(
        common_name=args.organism_name,
        taxon=args.taxon,
        ebi_id=args.ebi_id
    ))

    refseqs = {}
    for record in SeqIO.parse(args.fasta, "fasta"):
        refseq = post(token, args.url + 'refseq/', dict(  # noqa
            name=record.id,
            length=len(record.seq),
            organism=organism['id'],
        ))
        refseqs[record.id] = refseq

    for rec in GFF.parse(args.gff3):
        # rs = RefSeq.objects.get(name=rec.id, organism=organism)
        rs = refseqs[rec.id]
        for feat in rec.features:
            if feat.type != 'gene':
                continue

            post(token, args.url + 'genes/', {
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

if __name__ == '__main__':
    main()
