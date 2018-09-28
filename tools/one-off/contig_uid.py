#!/usr/bin/env python
import argparse
import hashlib
import multiprocessing.pool
from Bio import SeqIO
from Bio.Seq import reverse_complement

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def process_data(contig):
    data = []
    for i in range(len(contig)):
        # Slice beginning of sequence and append to the back
        rotation = "".join([contig[i:len(contig)], contig[0:i]])
        # Hash data and add to list
        data.append(hashlib.md5(rotation).hexdigest())
    return data


def generateUID(contig):
    #Take the contig and generate a list of each possible rotation of the sequence, and revcoms of each rotation
    #input is Bio.Seq object

    # Build a multiprocessing pool, one thread for forward, one for reverse
    pool = multiprocessing.pool.ThreadPool(processes=2)
    # And generate hashes of all rotations
    hashes = pool.map(process_data, [contig, reverse_complement(contig)], chunksize=1)
    pool.close()

    # Flatten this nested list of lits [[a, b], [c, d]] into [a, b, c, d]
    hashes = [item for sublist in hashes for item in sublist]

    # Sort
    hashes.sort()

    #Combine all hashes into a single digest to generate the UID
    uid = hashlib.md5()
    for h in hashes:
        uid.update(h)

    return uid.hexdigest()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a sequence UID for a contig')
    parser.add_argument('fa_files', type=argparse.FileType("r"), nargs="+", help="fasta files")

    args = parser.parse_args()
    for fa_file in args.fa_files:
        for record in SeqIO.parse(fa_file, 'fasta'):
            print("".join([record.id, "\t", generateUID(str(record.seq))]))
