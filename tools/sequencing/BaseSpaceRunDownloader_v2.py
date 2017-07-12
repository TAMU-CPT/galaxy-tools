#!/usr/bin/env python
from urllib2 import Request, urlopen, URLError
from os.path import expanduser
import re
import subprocess
import hashlib
import json
import sys
import os
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

CWD = os.getcwd()

CHUNK_SIZE = 1024 * 1024
RUN_ID = sys.argv[2]
API_BASE = 'https://api.basespace.illumina.com/'
AccessToken = sys.argv[1]
ACCESS_TOKEN = '?access_token=%s' % AccessToken
SampleListFilter = None
SampleListFilterEnabled = False
PARAMETERS = {'Limit': 1024, 'SortDir': 'Desc', 'SortBy': 'DateCreated'}

if len(sys.argv) > 3:
    with open(sys.argv[3], 'r') as handle:
        # Grab just the first column, expects at least line separated IDs and tab separated if there are additional cols
        SampleListFilter = [x.strip().split('\t')[0] for x in handle.readlines()]
        SampleListFilterEnabled = True


def restrequest(rawrequest):
    rawrequest = API_BASE + rawrequest + ACCESS_TOKEN
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

def restquery(rawrequest, parameters):
    # Parameters is a dict example: {'Limit': 50, 'SortDir': 'Asc'}
    # can be any of https://developer.basespace.illumina.com/docs/content/documentation/rest-api/api-reference#ResourceCollectionRequests
    rawrequest = API_BASE + rawrequest + ACCESS_TOKEN
    for query, value in parameters.items():
        rawrequest+= '&' + query + '=' + str(value)

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


# Provides access to a resource query that allows for limits and sorting parameters
json_obj = restquery('v1pre3/runs/%s/properties/Output.Samples/items' % RUN_ID, PARAMETERS)

for item in json_obj['Response']['Items']:
    if SampleListFilterEnabled and item['Content']['Id'] not in SampleListFilter:
        continue

    # Makes a request for the Sample resource
    sample_json_obj = restrequest(item['Content']['Href'])
    # Makes a request for the Sample's Files (Location of FASTQ file information)
    files_obj = restrequest(sample_json_obj['Response']['HrefFiles'])

    for file in files_obj['Response']['Items']:
        potential_filename = file['Name']
        safe_filename = re.sub('[^A-Za-z0-9._-]', '', potential_filename)
        safe_filename = os.path.join(CWD, 'output', safe_filename)
        fastq_url = file['HrefContent']
        log.info("Downloading %s", safe_filename)

        response = urlopen(API_BASE + fastq_url + ACCESS_TOKEN, timeout=1200)

        content_length = int(response.info().getheader('Content-Length'))
        expected_chunks = content_length / CHUNK_SIZE

        if os.path.exists(safe_filename) and os.path.getsize(safe_filename) == content_length:
            log.info("%s already downloaded and is correct size" % safe_filename)
            continue
        elif os.path.exists(safe_filename):
            log.warning("%s exists but may be a partial download / different file" % safe_filename)
            sys.exit()

        with open(safe_filename, 'wb') as handle:
            i = 0
            while True:
                i += 1
                chunk = response.read(CHUNK_SIZE)
                if not chunk:
                    break
                handle.write(chunk)
        log.info("Unzipping %s", safe_filename)
        subprocess.check_call(['gunzip', safe_filename])

print(json.dumps(json_obj, indent=2))
