#!/usr/bin/env python
import json
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('json', type=file, help='JSON Data')

    args = parser.parse_args()


    # https://cpt.tamu.edu/apollo/annotator/loadLink?loc=NC_005880:0..148317&organism=326&tracks=
    data = json.load(args.json)
    if len(data) > 1:
        raise Exception("More than one organism listed. TODO. Contact esr@tamu.edu")

    HTML_TPL = """
<html>
    <head>
        <title>Embedded Apollo Access</title>
        <style type="text/css">
            body {{
                    margin: 0;
            }}
            iframe {{
                border: 0;
                width: 100%;
                height: 100%
            }}
        </style>
    </head>
    <body>
         <iframe src="{base_url}/annotator/loadLink?loc={chrom}&organism={orgId}&tracklist=1"></iframe>
    </body>
</html>
    """

    print HTML_TPL.format(base_url=args.apollo, chrom="", orgId=data[0]['id'])
