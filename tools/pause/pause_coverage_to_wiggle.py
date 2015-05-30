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


def coverage_data(bam_file):
    data_f = ""
    data_r = ""
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
            map_f = {}
            map_r = {}

            for col in work_bam.fetch(chrom, start, end):
                #print " ".join(map(str, [col.qstart, col.qend, col.rlen, col.aend, col.alen, col.pos]))
                #   qstart   qend   rlen   aend    alen   pos
                #   0        145    145    13537   143    13394
                # reverse strand
                # start is  13395
                # end is 13537

                try:
                    for i in range(col.pos, col.aend):
                        if col.is_reverse:
                            if i not in map_r:
                                map_r[i] = 0
                            map_r[i] += 1
                        else:
                            if i not in map_f:
                                map_f[i] = 0
                            map_f[i] += 1
                except:
                    # HMmmmmmm.
                    pass
            # Write to file
            data_f += gen_header(bam_file.name, 'f')
            data_f += "variableStep chrom=%s\n" % chrom
            for i in range(start + 1, end + 1):
                if i in map_f:
                    data_f += "%s %.1f\n" % (i, map_f[i])
                else:
                    data_f += "%s 0.0\n" % i
            data_r += gen_header(bam_file.name, 'r')
            data_r += "variableStep chrom=%s\n" % chrom
            for i in range(start + 1, end + 1):
                if i in map_r:
                    data_r += "%s %.1f\n" % (i, map_r[i])
                else:
                    data_r += "%s 0.0\n" % i
    return (data_f, data_r)


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
                    'default': 'wig.coverage.f',
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
                    'default': 'wig.coverage.r',
                    'data_format': 'text/plain',
                    'default_format': 'TXT',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.pause2.coverage_to_wiggle',
            'appname': 'PAUSE2 BAM to coverage Wiggle',
            'appvers': '0.1',
            'appdesc': 'create wiggle file from coverage information',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    (data_f, data_r) = coverage_data(options['bam_file'])

    from galaxygetopt.outputfiles import OutputFiles
    off = OutputFiles(name='wig_f', GGO=opts)
    off.CRR(data=data_f)
    ofr = OutputFiles(name='wig_r', GGO=opts)
    ofr.CRR(data=data_r)
