#!/usr/bin/env python
# -*- coding: utf-8 -*-
import json
import sys
import requests
import argparse


def auth(creds, url):
    data = json.load(creds)["cacao"]
    r = requests.post(url + "api-token-auth/", data=data)
    return "JWT " + r.json()["token"]


def get(token, url, params=None):
    q = requests.get(url, headers={"Authorization": token}, params=params).json()
    return q


def post(token, url, data):
    q = requests.post(url, data=data, headers={"Authorization": token}).json()
    return q


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("url", help="Backend URL. Include trailing slash.")
    parser.add_argument(
        "creds", type=argparse.FileType("r"), help="json file with username/password"
    )

    parser.add_argument("--taxon", type=str, default="NA")
    parser.add_argument("--ebi_id", type=str, default="NA")
    parser.add_argument("organism_file", type=argparse.FileType("r"))

    args = parser.parse_args()

    token = auth(args.creds, args.url)

    organism_name = json.load(args.organism_file)[0]
    if "commonName" in organism_name:
        organism_name = organism_name["commonName"]
    elif "common_name" in organism_name:
        organism_name = organism_name["common_name"]
    else:
        raise Exception("Bad input organism data")

    get(token, args.url + "organisms/", dict(common_name=organism_name))["results"][0]

    refseqs = get(
        token, args.url + "refseq/", dict(organism__common_name=organism_name)
    )
    # TODO: multiple?
    refseq = refseqs["results"][0]

    gaf_data = []
    page = 1
    while True:
        gafs = get(
            token,
            args.url + "gafs/",
            dict(
                gene__refseq_id=refseq["id"],
                # gene__refseq_id='df8baca3-3af4-4d72-8045-6de08072be77',
                page=page,
            ),
        )
        if "detail" in gafs and gafs["detail"] == "Invalid page.":
            break

        for result in gafs["results"]:
            gaf_data.append(result)
        page += 1

    def get_go_term(id):
        if id.startswith("CPT:") or id.startswith("GO:"):
            try:
                r = requests.get("https://cpt.tamu.edu/onto_api/%s.json" % id)
                return r.json()["name"]
            except:
                return id
        else:
            return id

    fields = [
        "owner",
        "id",
        "db",
        "gene",
        "go_id",
        "db_reference",
        "evidence_code",
        "with_or_from",
        "aspect",
        "date",
        "assigned_by",
        "annotation_extension",
        "notes",
    ]

    def fix_field(gaf, f):
        if f == "owner":
            return gaf[f]["username"]
        elif f == "gene":
            return gaf[f]["id"]
        else:
            return gaf[f]

    header = "# " + ("\t".join(fields)) + "\tGO Term"
    print(header)
    import codecs

    UTF8Writer = codecs.getwriter("utf8")
    sys.stdout = UTF8Writer(sys.stdout)

    reload(sys)
    sys.setdefaultencoding("utf8")

    for gaf in gaf_data:
        fixed_fields = [fix_field(gaf, f).encode("utf-8") for f in fields]
        fixed_fields.append(get_go_term(gaf["go_id"]))
        ff = "\t".join(fixed_fields)
        print(unicode(ff).encode("utf-8"))


if __name__ == "__main__":
    main()
