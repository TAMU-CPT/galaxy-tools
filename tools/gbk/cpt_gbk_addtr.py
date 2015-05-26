#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import sys
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name='add_tr')
import copy

__doc__ = """
TR Addition
===========

This tools adds TR regions to a Genbank file.

The region you supply duplicated and added to the end of the genome. The
features from that region are likewise duplicated, and "\_rep" added to
their locus tag to indicate repeated features. If the specified TR
region ends inside a feature, it is not duplicated. This tool does
**not** renumber.

Requirements
------------

Your genomes MUST be opened such that the start of the terminal repeat
is the first base! This is unfortunately non-negotiable. If the
repeat\_region is in the middle of the genome, it will be duplicated
as-is, and appended to the end. This feature can be abused, but you
**must** be aware of it. See the section `"Feature
Abuse" <#feature-abuse>`__

You may specify repeats in one of the following ways:

-  Add a repeat\_region feature to the start of your genome in artemis
-  Specify a repeat end location

Either of these methods is sufficient to add terminal repeats.

Feature Abuse
-------------

Sometimes you may find yourself in the following situation:

-  You have a genome with a ``repeat_region``
-  During sequencing, this repeat region was not fully collapsed by the
   sequencer
-  The full repeat is at the start
-  A portion of the front of the repeat is also at the end (maybe this
   was done to improve annotation?)

In this event, you need to do the following:

-  Calculate how much of the repeat region was duplicated
-  Create the ``repeat_region`` feature at the start of the genome
-  Move the feature up by the number of duplicated bases

Then, when the tool is run, this subset of the (true) repeat region will
be duplicated and features added to the end. After that, you will need
to fix the location of both repeat regions to reflect the true repeat
region
"""


def add_tr(gbk_file=None, end=None, **kwd):
    from Bio import SeqIO
    record = list(SeqIO.parse(gbk_file, "genbank"))[0]
    cut_start = 0
    cut_end = 0

    # Find region of interest
    if end is None:
        repeat_regions = [x for x in record.features if x.type == 'repeat_region']
        if len(repeat_regions) == 1:
            cut_start = repeat_regions[0].location.start
            cut_end = repeat_regions[0].location.end
        elif len(repeat_regions) == 0:
            log.error("No repeats and end not specified")
            sys.exit(1)
        else:
            cut_start = repeat_regions[0].location.start
            cut_end = repeat_regions[0].location.end
            log.warning("Multiple repeat_regions found, using first")
    else:
        cut_start = 0
        cut_end = int(end)

    # Clone features in region of interest for duplication, this will actually
    # grab the repeat_region as well, conveniently
    clonefeats = [copy.deepcopy(x) for x in record.features if
                    cut_start <= x.location.start <= cut_end and
                    cut_start <= x.location.end <= cut_end]

    # For each cloned feature update the location
    for feat in clonefeats:
        feat.location += len(record.seq)
        if 'locus_tag' in feat.qualifiers:
            feat.qualifiers['locus_tag'][0] += '_rep'

        if 'note' in feat.qualifiers:
            feat.qualifiers['note'].append('CPT_TR dupcliated feature')
        else:
            feat.qualifiers['note'] = ['CPT_TR dupcliated feature']

    # Append extra sequence
    record.seq = record.seq + record.seq[cut_start:cut_end]
    # Add extended feature set
    record.features.extend(clonefeats)
    # Return
    return [record]


def passthrough(cb):
    opts = GGO(
        options=[
            ['gbk_file', 'Input Genome',
             {'required': True, 'validate': 'File/Input'}],
            ['end', 'End of Terminal Repeat Region. READ THE DIRECTIONS BELOW',
             {'validate': 'String', 'required': False}],
        ],
        outputs=[
            [
                'gbk',
                'Genbank with TR',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'gbk',
                    'data_format': 'genomic/annotated',
                    'default_format': 'Genbank',
                }
            ],
        ],
        defaults={
            'appid': 'edu.tamu.cpt.gbk.add_tr',
            'appname': 'Add Terminal Repeats',
            'appvers': '1.0',
            'appdesc': 'to a Genbank file with one or more genomes',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    return (opts, options, cb(**options))


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    (opts, options, result) = passthrough(add_tr)

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='gbk', GGO=opts)
    of.CRR(data=result)
