#!/usr/bin/env python
import argparse
import copy
import logging
import re
import sys
from Bio import SearchIO
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name='blastxml2gff3')

