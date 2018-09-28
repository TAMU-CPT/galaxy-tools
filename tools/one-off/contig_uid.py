#!/usr/bin/env python
import argparse
import hashlib
import multiprocessing.pool
from Bio import SeqIO
from Bio.Seq import reverse_complement

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def process_data(sequence):
    return hashlib.md5(sequence).hexdigest()


def generateUID(contig):
    #Take the contig and generate a list of each possible rotation of the sequence, and revcoms of each rotation
    #input is Bio.Seq object
    data = []
    #generate rotations and revcom rotations and hash

    for i in range(len(contig)):
        #Slice beginning of sequence and append to the back
        rotation = "".join([contig[i:len(contig)], contig[0:i]])
        #Revcom the rotated contig
        revcom = reverse_complement(rotation)
        #Hash both and add to hashes list
        data.append(rotation)
        data.append(revcom)

    pool = multiprocessing.pool.ThreadPool(processes=8)
    hashes = pool.map(process_data, data, chunksize=1)
    pool.close()

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
