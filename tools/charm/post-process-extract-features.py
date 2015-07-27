import fileinput
import sys

for line in fileinput.input():
    if line.startswith('##'):
        if 'gff-version' in line:
            print line,
        if '##FASTA' in line.upper():
            sys.exit(0)
    else:
        if '\tchromosome\t' in line or '\tannotation\t' in line:
            continue
        else:
            if '\tCDS\t' in line:
                line = line.replace('\tCDS\t', '\tpolypeptide\t')
            print line,
