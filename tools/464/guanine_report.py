#!/usr/bin/env python
import sys
import requests
import argparse
import json


def auth(creds, url):
    r = requests.post(url + "api-token-auth/", data=json.load(creds)["guanine"])
    return "JWT " + r.json()["token"]


def post_result(
    student_id, points_earned, points_possible, token, url, assessment_id, notes
):
    headers = {"Authorization": token}
    values = {
        "student": student_id,
        "assessment": assessment_id,
        "points_earned": points_earned,
        "points_possible": points_possible,
        "notes": notes,
    }
    r = requests.post(url + "results/", data=values, headers=headers)
    return r


def student_id(email, url, token):
    headers = {"Authorization": token}
    email = email.replace("@", "%40")
    student_url = url + "students/?email=" + email
    r = requests.get(student_url, headers=headers)
    try:
        return r.json()["results"][0]["id"]
    except:
        sys.stderr.write("Unknown student\n")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="post an assessment result")
    parser.add_argument("guanine_url", help="GUANINE Backend URL")
    parser.add_argument(
        "creds", type=argparse.FileType("r"), help="json file with username/password"
    )
    parser.add_argument("assessment_id", help="GUANINE Assessment ID")
    parser.add_argument("student_email", help="email of the student receiving result")
    parser.add_argument(
        "points_earned", type=int, help="how many points the student earned"
    )
    parser.add_argument(
        "points_possible",
        type=int,
        help="how many points were possible in the assessment",
    )
    parser.add_argument(
        "--notes", type=argparse.FileType("r"), help="Notes / metadata file"
    )
    args = parser.parse_args()

    token = auth(args.creds, args.guanine_url)
    sid = student_id(args.student_email, args.guanine_url, token)
    r = post_result(
        sid,
        args.points_earned,
        args.points_possible,
        token,
        args.guanine_url,
        args.assessment_id,
        args.notes.read() if args.notes else "",
    )
    if r.status_code in (200, 201):
        print("Success")
    else:
        sys.stderr.write("Failure: %s\n" % r.text)
        sys.exit(1)


if __name__ == "__main__":
    main()
