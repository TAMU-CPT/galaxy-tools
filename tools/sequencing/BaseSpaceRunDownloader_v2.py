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
import argparse

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

CWD = os.getcwd()
CHUNK_SIZE = 1024 * 1024
API_BASE = "https://api.basespace.illumina.com/"


def restrequest(rawrequest, access_token):

    ACCESS_TOKEN = "?access_token=%s" % access_token
    rawrequest = API_BASE + rawrequest + ACCESS_TOKEN

    log.info("Req: " + rawrequest.replace(access_token, "*" * len(access_token)))

    req_url_hash = hashlib.md5(rawrequest).hexdigest()
    cache_path_basedir = os.path.join(expanduser("~"), ".cache", "basespace")
    if not os.path.exists(cache_path_basedir):
        os.makedirs(cache_path_basedir)

    cache_path = os.path.join(cache_path_basedir, req_url_hash)
    if os.path.exists(cache_path):
        with open(cache_path, "r") as handle:
            json_obj = json.load(handle)
            return json_obj

    log.debug("Cache miss")
    request = Request(rawrequest)

    try:
        response = urlopen(request)
        json_string = response.read()
        json_obj = json.loads(json_string)
    except URLError as e:
        print("Got an error code:", e)
        sys.exit()

    with open(cache_path, "w") as handle:
        json.dump(json_obj, handle)

    return json_obj


def restquery(rawrequest, access_token, parameters):
    # Parameters is a dict example: {'Limit': 50, 'SortDir': 'Asc'}
    # can be any of https://developer.basespace.illumina.com/docs/content/documentation/rest-api/api-reference#ResourceCollectionRequests
    ACCESS_TOKEN = "?access_token=%s" % access_token
    rawrequest = API_BASE + rawrequest + ACCESS_TOKEN
    for query, value in parameters.items():
        rawrequest += "&" + query + "=" + str(value)

    log.info("Req: " + rawrequest.replace(access_token, "*" * len(access_token)))

    req_url_hash = hashlib.md5(rawrequest).hexdigest()
    cache_path_basedir = os.path.join(expanduser("~"), ".cache", "basespace")
    if not os.path.exists(cache_path_basedir):
        os.makedirs(cache_path_basedir)

    cache_path = os.path.join(cache_path_basedir, req_url_hash)
    if os.path.exists(cache_path):
        with open(cache_path, "r") as handle:
            json_obj = json.load(handle)
            return json_obj

    log.debug("Cache miss")
    request = Request(rawrequest)

    try:
        response = urlopen(request)
        json_string = response.read()
        json_obj = json.loads(json_string)
    except URLError as e:
        print("Got an error code:", e)
        sys.exit()

    with open(cache_path, "w") as handle:
        json.dump(json_obj, handle)

    return json_obj


def rundownload(
    access_token, runid, ids=None, limit=1024, offset=0, sortby="Id", sortdir="Asc"
):
    # Provides access to a resource query that allows for limits and sorting parameters
    parameters = {
        "Limit": limit,
        "Offset": offset,
        "SortBy": sortby,
        "SortDir": sortdir,
    }

    if ids:
        # Grab just the first column, expects at least line separated IDs and tab separated if there are additional cols
        SampleListFilter = [x.strip().split("\t")[0] for x in ids.readlines()]
        SampleListFilterEnabled = True

    json_obj = restquery(
        "v1pre3/runs/%s/properties/Output.Samples/items" % runid,
        access_token,
        parameters,
    )
    responses = [json_obj]

    if json_obj["Response"]["TotalCount"] > limit:
        parameters["Offset"] = limit
        json_obj2 = restquery(
            "v1pre3/runs/%s/properties/Output.Samples/items" % runid,
            access_token,
            parameters,
        )
        responses.append(json_obj2)

    ACCESS_TOKEN = "?access_token=%s" % access_token

    # test_count = 0
    for response in responses:
        for item in response["Response"]["Items"]:
            if SampleListFilterEnabled:
                if item["Content"]["Id"] not in SampleListFilter:
                    continue

            # Makes a request for the Sample resource
            sample_json_obj = restrequest(item["Content"]["Href"], access_token)
            # Makes a request for the Sample's Files (Location of FASTQ file information)
            files_obj = restrequest(
                sample_json_obj["Response"]["HrefFiles"], access_token
            )
            for file in files_obj["Response"]["Items"]:
                potential_filename = file["Name"]
                safe_filename = re.sub("[^A-Za-z0-9._-]", "", potential_filename)
                safe_filename = os.path.join(CWD, "output", safe_filename)
                fastq_url = file["HrefContent"]
                log.info("Downloading %s", safe_filename)

                response = urlopen(API_BASE + fastq_url + ACCESS_TOKEN, timeout=1200)

                content_length = int(response.info().getheader("Content-Length"))
                # expected_chunks = content_length / CHUNK_SIZE
                # expected_chunks meant to be used to check download but not implemented. Removing for pyflakes

                if (
                    os.path.exists(safe_filename)
                    and os.path.getsize(safe_filename) == content_length
                ):
                    log.info(
                        "%s already downloaded and is correct size" % safe_filename
                    )
                    continue
                elif os.path.exists(safe_filename):
                    log.warning(
                        "%s exists but may be a partial download / different file"
                        % safe_filename
                    )
                    sys.exit()

                with open(safe_filename, "wb") as handle:
                    i = 0
                    while True:
                        i += 1
                        chunk = response.read(CHUNK_SIZE)
                        if not chunk:
                            break
                        handle.write(chunk)
                log.info("Unzipping %s", safe_filename)
                subprocess.check_call(["gunzip", safe_filename])
            # Limit to 4 samples for testing
            # test_count += 1
            # if test_count > 4:
            #    break

    print(json.dumps(responses, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download all fastq files from BaseSpace for a given run"
    )
    parser.add_argument("access_token", type=str, help="BaseSpace API access Token")
    parser.add_argument("runid", type=str, help="Run ID")
    # Tabular list of sample IDs (Illumina IDs) to retrieve from Run. If left out, defaults to collect ALL samples in run
    parser.add_argument("--ids", type=argparse.FileType("r"), help="List of Sample IDs")

    # Exposes Collection Query arguments if necessary
    parser.add_argument(
        "--limit",
        nargs="?",
        type=int,
        help="Number of samples to retrieve",
        default=1024,
    )
    parser.add_argument(
        "--offset",
        nargs="?",
        type=int,
        help="Offset to begin retrieving samples",
        default=0,
    )
    parser.add_argument(
        "--sortby",
        nargs="?",
        type=str,
        help="Sort by Id or by DateCreated",
        default="Id",
    )
    parser.add_argument(
        "--sortdir",
        nargs="?",
        type=str,
        help="Sort (Asc)ending or (Desc)ending",
        default="Asc",
    )

    args = parser.parse_args()
    rundownload(**vars(args))
