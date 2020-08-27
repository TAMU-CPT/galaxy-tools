#!/usr/bin/env python
from __future__ import print_function

import argparse
import json


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('json', type=argparse.FileType("r"), help='JSON Data')
    parser.add_argument('external_apollo_url')
    parser.add_argument('--tracklist', action='store_true')

    args = parser.parse_args()

    # https://fqdn/apollo/annotator/loadLink?loc=NC_005880:0..148317&organism=326
    data = json.load(args.json)

    # This is base64 encoded to get past the toolshed's filters.
    HTML_TPL = """
        <html>
            <head>
                <title>Embedded Apollo Access</title>
                <style type="text/css">body {{margin: 0;}} iframe {{border: 0;width: 100%;height: 100%}}</style>
            </head>
            <body>
                <iframe src="{base_url}/annotator/loadLink?loc={chrom}&organism={orgId}{tracklist}"></iframe>
            </body>
        </html>
    """

    print(HTML_TPL.format(base_url=args.external_apollo_url, chrom="", orgId=data[0]['id'], tracklist='&tracklist=1' if args.tracklist else ''))
