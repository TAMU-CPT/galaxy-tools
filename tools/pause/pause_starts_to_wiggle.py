#!/usr/bin/env python
"""Extract read start from BAM files to Wig format for PAUSE.

Usage:
    bam_to_wiggle.py <BAM file>

"""
import os
from contextlib import contextmanager
import pysam
from galaxygetopt.ggo import GalaxyGetOpt as GGO


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


def start_data(bam_file):
    wig_f = ""
    wig_r = ""
    with indexed_bam(bam_file) as work_bam:
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
                #print " ".join(map(str, [col.qstart, col.qend, col.rlen, col.aend, col.alen, col.pos]))
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
            wig_f += gen_header(bam_file.name, 'f')
            wig_f += "variableStep chrom=%s\n" % chrom
            for i in range(start + 1, end + 1):
                if i in start_map_f:
                    wig_f += "%s %.1f\n" % (i, start_map_f[i])
                else:
                    wig_f += "%s 0.0\n" % i
            wig_r += gen_header(bam_file.name, 'r')
            wig_r += "variableStep chrom=%s\n" % chrom
            for i in range(start + 1, end + 1):
                if i in start_map_r:
                    wig_r += "%s %.1f\n" % (i, start_map_r[i])
                else:
                    wig_r += "%s 0.0\n" % i
    return (wig_f, wig_r)


if __name__ == "__main__":
    opts = GGO(
        options=[
            ['bam_file', 'Bam File',
             {'required': True, 'validate': 'File/Input'}],
        ],
        outputs=[
            [
                'wig_f',
                '+ strand wig data',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'wig.starts.f',
                    'data_format': 'text/plain',
                    'default_format': 'TXT',
                }
            ],
            [
                'wig_r',
                '- strand wig data',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'wig.starts.r',
                    'data_format': 'text/plain',
                    'default_format': 'TXT',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.pause2.starts_to_wiggle',
            'appname': 'PAUSE2 BAM to Starts Wiggle',
            'appvers': '0.1',
            'appdesc': 'create wiggle file from starts information',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    (data_f, data_r) = start_data(options['bam_file'])

    from galaxygetopt.outputfiles import OutputFiles
    off = OutputFiles(name='wig_f', GGO=opts)
    off.CRR(data=data_f)
    ofr = OutputFiles(name='wig_r', GGO=opts)
    ofr.CRR(data=data_r)
