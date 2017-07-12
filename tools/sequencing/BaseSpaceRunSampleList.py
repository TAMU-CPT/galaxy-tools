#!/usr/bin/env python
from urllib2 import Request, urlopen, URLError
from os.path import expanduser
import hashlib
import json
import sys
import os
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

CHUNK_SIZE = 128 * 1024
RUN_ID = sys.argv[2]
API_BASE = 'https://api.basespace.illumina.com/'
AccessToken = sys.argv[1]
ACCESS_TOKEN = '?access_token=%s' % AccessToken


def restrequest(rawrequest):
    rawrequest = API_BASE + rawrequest + ACCESS_TOKEN
    log.debug('Req: ' + rawrequest.replace(AccessToken, '*' * len(AccessToken)))

    req_url_hash = hashlib.md5(rawrequest).hexdigest()
    cache_path_basedir = os.path.join(expanduser("~"), '.cache', 'basespace')
    if not os.path.exists(cache_path_basedir):
        os.makedirs(cache_path_basedir)

    cache_path = os.path.join(cache_path_basedir, req_url_hash)
    if os.path.exists(cache_path):
        with open(cache_path, 'r') as handle:
            json_obj = json.load(handle)
            return json_obj

    log.debug("Cache miss")
    request = Request(rawrequest)

    try:
        response = urlopen(request)
        json_string = response.read()
        json_obj = json.loads(json_string)
    except URLError as e:
        print('Got an error code:', e)
        sys.exit()

    with open(cache_path, 'w') as handle:
        json.dump(json_obj, handle)

    return json_obj


def restquery(rawrequest, parameters):
    # Parameters is a dict example: {'Limit': 50, 'SortDir': 'Asc'}
    # can be any of https://developer.basespace.illumina.com/docs/content/documentation/rest-api/api-reference#ResourceCollectionRequests
    rawrequest = API_BASE + rawrequest + ACCESS_TOKEN
    for query, value in parameters.items():
        rawrequest += '&' + query + '=' + str(value)

    log.info('Req: ' + rawrequest.replace(AccessToken, '*' * len(AccessToken)))

    req_url_hash = hashlib.md5(rawrequest).hexdigest()
    cache_path_basedir = os.path.join(expanduser("~"), '.cache', 'basespace')
    if not os.path.exists(cache_path_basedir):
        os.makedirs(cache_path_basedir)

    cache_path = os.path.join(cache_path_basedir, req_url_hash)
    if os.path.exists(cache_path):
        with open(cache_path, 'r') as handle:
            json_obj = json.load(handle)
            return json_obj

    log.debug("Cache miss")
    request = Request(rawrequest)

    try:
        response = urlopen(request)
        json_string = response.read()
        json_obj = json.loads(json_string)
    except URLError as e:
        print('Got an error code:', e)
        sys.exit()

    with open(cache_path, 'w') as handle:
        json.dump(json_obj, handle)

    return json_obj


json_obj = restquery('v1pre3/runs/%s/properties/Output.Samples/items' % RUN_ID, {'Limit': 1024})
responses = [json_obj]

if json_obj['Response']['TotalCount'] > 1024:
    newoffset = json_obj['Response']['TotalCount'] - 1024
    json_obj2 = restquery('v1pre3/runs/%s/properties/Output.Samples/items' % RUN_ID, {'Limit': 1024, 'Offset': newoffset})
    responses.append(json_obj2)

print(json.dumps(obj) for obj in responses)
