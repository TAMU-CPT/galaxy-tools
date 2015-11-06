#!/usr/bin/env python
import argparse
import requests


def addOrganism(apollo, username, password, cn, jbrowse, blatdb=None,
                genus=None, species=None):

    payload = {
        'username': username,
        'password': password,

        'organism': {
            'commonName': cn,
            'directory': jbrowse,
            'blatdb': blatdb,
            'genus': genus,
            'species': species,
        }
    }
    r = requests.post(apollo + '/organism/addOrganism', data=payload)
    print(r.text)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Username')
    parser.add_argument('password', help='WA Password')

    parser.add_argument('cn', help='Organism Common Name')
    parser.add_argument('jbrowse', help='JBrowse Data Directory')
    parser.add_argument('--blatdb', help='BlatDB Directory')
    parser.add_argument('--genus', help='Organism Genus')
    parser.add_argument('--species', help='Organism Species')

    args = parser.parse_args()
    addOrganism(**vars(args))
