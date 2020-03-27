#!/usr/bin/env python
import json
import argparse
from webapollo import WAAuth, WebApolloInstance, AssertUser, accessible_organisms

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="List all organisms available in an Apollo instance"
    )
    WAAuth(parser)
    parser.add_argument("email", help="User Email")
    args = parser.parse_args()

    wa = WebApolloInstance(args.apollo, args.username, args.password)

    gx_user = AssertUser(wa.users.loadUsers(email=args.email))
    all_orgs = wa.organisms.findAllOrganisms()

    orgs = accessible_organisms(gx_user, all_orgs)

    cleanedOrgs = []
    for organism in all_orgs:
        org = {
            "name": organism["commonName"],
            "id": organism["id"],
            "annotations": organism["annotationCount"],
            "sequences": organism["sequences"],
        }
        cleanedOrgs.append(org)

    print(json.dumps(cleanedOrgs, indent=2))
