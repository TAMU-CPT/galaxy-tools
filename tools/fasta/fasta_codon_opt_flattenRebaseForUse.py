#!/usr/bin/env python
import sys
import yaml

d = yaml.load(sys.stdin)
q = {}

for x in d:
    q[x] = [w[3:] for w in d[x]["recognition_sequence"]]

print(yaml.dump(q))
