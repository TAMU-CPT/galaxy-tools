#!/usr/bin/env python
import argparse
import logging
import sys

logging.basicConfig(level=logging.INFO)


def extract_clust(cdhit):
    cluster = []
    cluster_id = None
    for line in cdhit:
        if line.startswith(">"):
            yield cluster, cluster_id
            cluster = []
            cluster_id = line.strip().split(" ")[1]
        elif cluster_id is not None:
            protein_id = line[line.index(", >") + 3 : line.rindex("... ")]
            # Move representative to front.
            if "... *" in line:
                cluster = [protein_id] + cluster
            else:
                cluster.append(protein_id)

    yield cluster, cluster_id


def cdhit_process(cdhit):
    for (cluster_elements, cluster_id) in extract_clust(cdhit):
        if cluster_id is not None:
            sys.stdout.write(cluster_id + "\t" + str(len(cluster_elements)) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cdhit", type=argparse.FileType("r"))
    args = parser.parse_args()
    cdhit_process(**vars(args))
