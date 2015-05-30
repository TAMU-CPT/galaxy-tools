#!/usr/bin/env python
import os
import sys
import argparse
import subprocess
import tempfile
from contextlib import contextmanager, closing


def main(xmfa_file, window_size=100, relative_to=1, sequences=None):
    config = {"program": {"ucsc_bedGraphToBigWig": ["bedGraphToBigWig"]}}
    if outfile is None:
        outfile = "%s.bigwig" % os.path.splitext(bam_file)[-1]


    sizes = get_sizes(bam_file, config)
    print "Have %i references" % len(sizes)
    if not sizes:
        sys.stderr.write("Problem reading BAM header.\n")
        sys.exit(1)

    #Use a temp file to avoid any possiblity of not having write permission
    temp_handle = tempfile.NamedTemporaryFile(delete=False)
    temp_file = temp_handle.name
    with closing(temp_handle):
        print "Calculating coverage..."
        convert_to_graph(bam_file, split, config, temp_handle)
    try:
        print "Converting %i MB graph file to bigwig..." % (os.path.getsize(temp_file) // (1024 * 1024))
        #Can't pipe this as stdin due to converter design,
        #https://lists.soe.ucsc.edu/pipermail/genome/2011-March/025455.html
        convert_to_bigwig(temp_file, sizes, config, outfile)
    finally:
        if os.path.isfile(temp_file):
            os.remove(temp_file)
    print "Done"

@contextmanager
def indexed_bam(bam_file, config):
    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file)
    sam_reader = pysam.Samfile(bam_file, "rb")
    yield sam_reader
    sam_reader.close()

def get_sizes(bam_file, config):
    with indexed_bam(bam_file, config) as work_bam:
        sizes = zip(work_bam.references, work_bam.lengths)
    return sizes

def convert_to_graph(bam_file, split, config, out_handle):
    cl = config["program"]["bedtools_genomeCoverageBed"] + ["-ibam", bam_file, "-bg"]
    if split:
        cl.append("-split")
    subprocess.check_call(cl, stdout=out_handle)

def convert_to_bigwig(bedgraph_file, chr_sizes, config, bw_file):
    #This will be fine under Galaxy, but could use temp folder?
    size_file = "%s-sizes.txt" % (os.path.splitext(bw_file)[0])
    with open(size_file, "w") as out_handle:
        for chrom, size in chr_sizes:
            out_handle.write("%s\t%s\n" % (chrom, size))
    try:
        cl = config["program"]["ucsc_bedGraphToBigWig"] + [bedgraph_file, size_file, bw_file]
        subprocess.check_call(cl)
    finally:
        os.remove(size_file)
    return bw_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert XMFA alignments to gff3', prog='xmfa2gff3')
    parser.add_argument('xmfa_file', type=file, help='XMFA File')
    parser.add_argument('--window_size', type=int, help='Window size for analysis', default=1000)
    parser.add_argument('--relative_to', type=str, help='Index of the parent sequence in the MSA', default='1')
    parser.add_argument('--sequences', type=file, nargs='+',
                        help='Fasta files (in same order) passed to parent for reconstructing proper IDs')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()
    main(**vars(args))
