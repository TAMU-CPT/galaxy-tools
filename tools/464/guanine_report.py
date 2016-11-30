import requests
import argparse
import json


def auth(creds):
    url = "http://localhost:8000/api-token-auth/"
    r = requests.post(url, data=json.load(creds))
    return 'JWT ' + r.json()['token']


def post_result(student_id, points_earned, points_possible, token):
    url = "http://localhost:8000/results/"
    headers = {'Authorization': token}
    values = {'student': student_id,
              'assessment': 'd18cfd41-85e2-4e75-8333-22339e05edc2',
              'points_earned': points_earned,
              'points_possible': points_possible}
    r = requests.post(url, data=values, headers=headers)
    return r


def student_id(email):
    email = email.replace('@', '%40')
    url = "http://localhost:8000/students/?email=" + email
    r = requests.get(url)
    return r.json()['results'][0]['id']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='post an assessment result')
    parser.add_argument('creds', type=file, help='json file with username/password')
    parser.add_argument('student_email', help='email of the student receiving result')
    parser.add_argument('points_earned', type=int, help='how many points the student earned')
    parser.add_argument('points_possible', type=int, help='how many points were possible in the assessment')
    args = parser.parse_args()

    token = auth(args.creds)
    student_id = student_id(args.student_email)
    post_result(student_id, args.points_earned, args.points_possible, token)
