#!/usr/bin/env python
import argparse
import requests


def addAttribute(webapollo, username, password, organism_cn, track,
                 unique_name, attribute, value):

    payload = {
        'username': username,
        'password': password,
        'features': [
            {
                'non_reserved_properties': [
                    {
                        "tag": attribute,
                        "value": value,
                    }
                ],
                "uniquename": unique_name,
            }
        ],
        'track': track,
        'organism': organism_cn
    }
    r = requests.post(webapollo + '/annotationEditor/addAttribute', data=payload)
    print(r.text)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('webapollo', help='Complete WebApollo URL')
    parser.add_argument('username', help='WA Username')
    parser.add_argument('password', help='WA Password')
    parser.add_argument('organism_cn', help='Organism Common Name')
    parser.add_argument('track', help='Track Name')
    parser.add_argument('unique_name', help='Unique name for feature')
    parser.add_argument('attribute', help='Attribute name')
    parser.add_argument('value', help='Attribute Value')

    args = parser.parse_args()
    addAttribute(**vars(args))
