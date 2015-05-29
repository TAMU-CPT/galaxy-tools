#!/usr/bin/env python
"""Convert BAM files to BigWig file format in a specified region.

Original version copyright Brad Chapman with revisions from Peter Cock
and ideas from Lance Parsons

Usage:
    bam_to_bigwig.py <BAM file> [--outfile=<output file name>] [--split]

The --split argument is passed to bedtools genomecov

The script requires:
    bedGraphToBigWig from UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/)
"""
import os
import sys
import subprocess
import tempfile
from optparse import OptionParser
from contextlib import contextmanager, closing

import pysam

def main(bam_file, outfile=None, split=False):
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
    parser = OptionParser()
    parser.add_option("-o", "--outfile", dest="outfile")
    parser.add_option("-s", "--split", action="store_true", dest="split")
    (options, args) = parser.parse_args()
    if len(args) not in [1, 2]:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict(
        outfile=options.outfile,
        split=options.split)
    main(*args, **kwargs)
