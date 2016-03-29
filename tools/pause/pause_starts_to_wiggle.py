#!/usr/bin/env python
"""Extract read start from BAM files to Wig format for PAUSE.

Usage:
    bam_to_wiggle.py <BAM file>

"""
import os
import tempfile
from contextlib import contextmanager
import pysam
import subprocess
import argparse


@contextmanager
def indexed_bam(bam_file):
    if not os.path.exists(bam_file.name + ".bai"):
        pysam.index(bam_file.name)
    sam_reader = pysam.Samfile(bam_file.name, "rb")
    yield sam_reader
    sam_reader.close()


def gen_header(bam_file, suffix):
    track_name = "name=%s_%s" % (os.path.splitext(
        os.path.split(bam_file)[-1])[0], suffix)
    return "track type=wiggle_0 %s visibility=full\n" % track_name


def convert_to_bigwig(wig_file, chr_sizes, bw_file):
    # This will be fine under Galaxy, but could use temp folder?
    size_file = "%s-sizes.txt" % (os.path.splitext(bw_file)[0])
    with open(size_file, "w") as out_handle:
        for chrom, size in chr_sizes:
            out_handle.write("%s\t%s\n" % (chrom, size))
    try:
        cl = ["wigToBigWig", wig_file, size_file, bw_file]
        subprocess.check_call(cl)
    finally:
        os.unlink(size_file)
    return bw_file


def start_data(bam_file, starts_f=None, starts_r=None):
    with indexed_bam(bam_file) as work_bam:
        starts_f_wig = tempfile.NamedTemporaryFile(delete=False)
        starts_r_wig = tempfile.NamedTemporaryFile(delete=False)

        sizes = zip(work_bam.references, work_bam.lengths)
        regions = [(name, 0, length) for name, length in sizes]
        for chrom, start, end in regions:
            if end is None and chrom in work_bam.references:
                end = work_bam.lengths[work_bam.references.index(chrom)]
            assert end is not None, "Could not find %s in header" % chrom

            # Since the file is sorted, we could actually optimise this bit
            # out...currently fails cost benefit analysis so will wait until
            # memory issues are reported.
            start_map_f = {}
            start_map_r = {}

            for col in work_bam.fetch(chrom, start, end):
                # print " ".join(map(str, [col.qstart, col.qend, col.rlen, col.aend, col.alen, col.pos]))
                #   qstart   qend   rlen   aend    alen   pos
                #   0        145    145    13537   143    13394
                # reverse strand
                # start is  13395
                # end is 13537
                if col.is_reverse:
                    rstart = col.aend
                    if rstart in start_map_r:
                        start_map_r[rstart] += 1
                    else:
                        start_map_r[rstart] = 1
                else:
                    rstart = col.pos + 1
                    if rstart in start_map_f:
                        start_map_f[rstart] += 1
                    else:
                        start_map_f[rstart] = 1
            # Write to file
            starts_f_wig.write(gen_header(bam_file.name, 'f'))
            starts_f_wig.write("variableStep chrom=%s\n" % chrom)
            for i in range(start + 1, end + 1):
                if i in start_map_f:
                    starts_f_wig.write("%s %.1f\n" % (i, start_map_f[i]))
                else:
                    starts_f_wig.write("%s 0.0\n" % i)
            starts_r_wig.write(gen_header(bam_file.name, 'r'))
            starts_r_wig.write("variableStep chrom=%s\n" % chrom)
            for i in range(start + 1, end + 1):
                if i in start_map_r:
                    starts_r_wig.write("%s %.1f\n" % (i, start_map_r[i]))
                else:
                    starts_r_wig.write("%s 0.0\n" % i)

        starts_f_wig.close()
        starts_r_wig.close()

        try:
            convert_to_bigwig(starts_f_wig.name, sizes, starts_f.name)
            convert_to_bigwig(starts_r_wig.name, sizes, starts_r.name)
        finally:
            os.unlink(starts_f_wig.name)
            os.unlink(starts_r_wig.name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract starts from BAM as BigWig')
    parser.add_argument('bam_file', type=file, help='Bam file')
    parser.add_argument('--starts_f', type=argparse.FileType('wb'), default='starts.f.bw', help='Sense Starts File')
    parser.add_argument('--starts_r', type=argparse.FileType('wb'), default='starts.r.bw', help='Antisense Starts File')
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()

    start_data(**vars(args))
