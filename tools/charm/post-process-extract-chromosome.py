import fileinput
import sys

for line in fileinput.input():
    if line.startswith('##'):
        if 'gff-version' in line:
            print line,
        if '##FASTA' in line.upper():
            sys.exit(0)
    else:
        if '\tchromosome\t' in line:
            print line,
        elif '\tannotation\t' in line:
            continue
        else:
            sys.exit(0)
