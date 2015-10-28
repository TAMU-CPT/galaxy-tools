from Bio import SeqIO
import math
import sys
import logging
logging.basicConfig()
log = logging.getLogger()

data = {}
for record in SeqIO.parse(sys.argv[1], 'fasta'):
    data[record.id] = len(record)

for row in sys.stdin:
    rowdata = row.strip().split(' ')

    if rowdata[0] not in data:
        log.error("Unknown contig in bam not found in fasta")
        continue

    dlen = data[rowdata[0]]
    dsum = int(rowdata[1])
    dsumsq = int(rowdata[2])

    avg = dsum / dlen
    stdev = math.sqrt(dsumsq/dlen - math.pow((dsum/dlen), 2))

    print '%s\t%s\t%s' % (rowdata[0], avg, stdev)
