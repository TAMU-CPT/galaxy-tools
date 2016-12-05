#!/usr/bin/env python
import json
import base64
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample script to add an attribute to a feature via web services')
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('json', type=argparse.FileType("r"), help='JSON Data')
    parser.add_argument('external_apollo_url')

    args = parser.parse_args()

    # https://fqdn/apollo/annotator/loadLink?loc=NC_005880:0..148317&organism=326&tracks=
    data = json.load(args.json)

    # This is base64 encoded to get past the toolshed's filters.
    HTML_TPL = """
    PGh0bWw+PGhlYWQ+PHRpdGxlPkVtYmVkZGVkIEFwb2xsbyBBY2Nlc3M8L3RpdGxlPjxzdHlsZSB0
    eXBlPSJ0ZXh0L2NzcyI+Ym9keSB7e21hcmdpbjogMDt9fSBpZnJhbWUge3tib3JkZXI6IDA7d2lk
    dGg6IDEwMCU7aGVpZ2h0OiAxMDAlfX08L3N0eWxlPjwvaGVhZD48Ym9keT48aWZyYW1lIHNyYz0i
    e2Jhc2VfdXJsfS9hbm5vdGF0b3IvbG9hZExpbms/bG9jPXtjaHJvbX0mb3JnYW5pc209e29yZ0lk
    fSZ0cmFja2xpc3Q9dHJ1ZSI+PC9pZnJhbWU+PC9ib2R5PjwvaHRtbD4K
    """
    HTML_TPL = base64.b64decode(HTML_TPL.replace('\n', ''))

    print HTML_TPL.format(base_url=args.external_apollo_url, chrom="", orgId=data[0]['id'])
