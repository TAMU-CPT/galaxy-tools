#!/usr/bin/env python
import argparse
import logging
from graphviz import Digraph

logging.basicConfig(level=logging.INFO)


def xmfa_graph(xmfa_backbone):
    dot = Digraph(comment="XMFA")
    nodes = [{"id": "start", "text": "Alignment Start"}]
    edges = []
    last_id = "start"
    node_genome_last_seen = {}

    parsed = []
    for idx, line in enumerate(xmfa_backbone):
        if line.startswith("seq"):
            continue
        parsed = map(int, line.split("\t"))

        # Setup last-seen data
        if len(node_genome_last_seen.keys()) == 0:
            for i in range(len(parsed) / 2):
                node_genome_last_seen[i] = "start"
        last_id = "gen_%s" % idx

        text = ""
        for i in range(len(parsed) / 2):
            if parsed[2 * i] != 0:
                text += "Seq %s: %s-%s\n" % (i, parsed[2 * i], parsed[2 * i + 1])
                edges.append(
                    {"from": node_genome_last_seen[i], "to": last_id, "txt": "%s" % i}
                )
                node_genome_last_seen[i] = last_id
        nodes.append({"id": last_id, "text": text})

    nodes.append({"id": "end", "text": "Alignment End"})

    for i in range(len(parsed) / 2):
        edges.append({"from": node_genome_last_seen[i], "to": "end", "txt": "%s" % i})
        node_genome_last_seen[i] = "end"

    # Create dot plot
    for node in nodes:
        dot.node(node["id"], node["text"])

    for edge in edges:
        dot.edge(edge["from"], edge["to"], label=edge["txt"])

    dot.render("xmfa.gv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Graph of XMFA backbone")
    parser.add_argument(
        "xmfa_backbone", type=argparse.FileType("r"), help="XMFA Backbone"
    )

    args = parser.parse_args()
    xmfa_graph(**vars(args))
