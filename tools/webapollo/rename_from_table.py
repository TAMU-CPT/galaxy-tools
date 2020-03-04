#!/usr/bin/env python
import argparse
import sys
from webapollo import WebApolloInstance
from webapollo import WAAuth, OrgOrGuess, GuessOrg, AssertUser, retry
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    WAAuth(parser)
    parser.add_argument("email", help="User Email")
    parser.add_argument("table", help="Data Table")
    parser.add_argument(
        "--set", choices=["name", "qualifier"], help="Which attribute to set."
    )
    parser.add_argument("--qualifier_name", default="locus_tag")
    OrgOrGuess(parser)

    args = parser.parse_args()

    with open(args.table, "r") as f:
        # reads table into data dict ID: Value
        data = {}
        for line in f:
            (key, val) = line.split("\t")
            data[key] = val.strip()

    wa = WebApolloInstance(args.apollo, args.username, args.password)
    # User must have an account
    gx_user = AssertUser(wa.users.loadUsers(email=args.email))

    # Get organism
    org_cn = GuessOrg(args, wa)
    if isinstance(org_cn, list):
        org_cn = org_cn[0]

    # TODO: Check user perms on org.
    org = wa.organisms.findOrganismByCn(org_cn)

    # Call setSequence to tell apollo which organism we're working with
    wa.annotations.setSequence(org["commonName"], org["id"])
    # Then get a list of features.
    features = wa.annotations.getFeatures()
    # For each feature in the features
    sys.stdout.write("# Feature ID\tNew Name\n")
    for feature in sorted(features["features"], key=lambda x: x["location"]["fmin"]):
        new_name = None
        if "parent_id" in feature:
            if feature["parent_id"] in data:
                new_name = data[feature["parent_id"]]
            elif feature["uniquename"] in data:
                new_name = data[feature["uniquename"]]
            else:
                sys.stdout.write("%s\tNot Found\n" % feature["uniquename"])
                continue
        else:
            if feature["uniquename"] in data:
                new_name = data[feature["uniquename"]]
            else:
                sys.stdout.write("%s\tNot Found\n" % feature["uniquename"])
                continue
        if args.set == "name":

            def fn0():
                wa.annotations.setName(feature["parent_id"], new_name)

            def fn1():
                wa.annotations.setName(feature["uniquename"], new_name)

            retry(fn0)
            retry(fn1)
        else:

            def fn0():
                wa.annotations.addAttributes(
                    feature["parent_id"], {"locus_tag": [new_name]}
                )

            retry(fn0)
        sys.stdout.write("%s\t%s\n" % (feature["uniquename"], new_name))
