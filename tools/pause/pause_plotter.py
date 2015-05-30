#!/usr/bin/env python
"""PAUSE: Plotter
"""
import bx.wiggle
import numpy
import cpt_pause
from galaxygetopt.ggo import GalaxyGetOpt as GGO


def get_data(wig_handle, count):
    data = []
    for row in bx.wiggle.Reader(wig_handle):
        data.append((row[1], row[2]))

    reshaped = numpy.array(data, dtype=numpy.int)
    # Odd rows need to be fixed
    if count % 2 == 1:
        reshaped[:, 1] *= -1
    return reshaped


def main(coverage=None, starts=None, highlights=None):
    if coverage is None:
        coverage = []
    if starts is None:
        starts = []
    if highlights is None:
        highlights = []

    track_list = []
    count = 0

    # Coverage is handled separately and "just for looks"
    for wig_handle in coverage:
        # Fetch data
        reshaped = get_data(wig_handle, count)
        # Downsample
        reshaped = cpt_pause.Filter.downsample(reshaped, sampling_interval=10)
        # Append
        track_list.append(cpt_pause.Coverage(reshaped, opacity=0.5))
        count += 1

    # Starts are handled separately from coverage
    for wig_handle in starts:
        # Fetch data
        reshaped = get_data(wig_handle, count)
        # Pass filter, then remove repeated y values
        reshaped = cpt_pause.Filter.repeat_reduction(
            cpt_pause.Filter.minpass(reshaped, min_value=2))
        track_list.append(cpt_pause.Coverage(reshaped, line_color='blue'))
        count += 1

    for wig_handle in highlights:
        # Fetch data
        reshaped = get_data(wig_handle, count)
        track_list.append(cpt_pause.Highlight(reshaped))
        count += 1

    g = cpt_pause.Gfx(track_list)
    return g.plot()


if __name__ == "__main__":
    opts = GGO(
        options=[
            ['coverage', 'Coverage files',
             {'multiple': True, 'validate': 'File/Input'}],
            ['starts', 'Start files',
             {'multiple': True, 'validate': 'File/Input'}],
            ['highlights', 'Data Highlights',
             {'multiple': True, 'validate': 'File/Input'}],
        ],
        outputs=[
            [
                'pause_plot',
                'PAUSE Plot',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'pause_plot',
                    'data_format': 'text/plain',
                    'default_format': 'TXT',
                }
            ]
        ],
        defaults={
            'appid': 'edu.tamu.cpt.pause2.plotter',
            'appname': 'PAUSE2 Plotter',
            'appvers': '0.3',
            'appdesc': 'run PAUSE plotting',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    output = main(coverage=options['coverage'], starts=options['starts'],
                  highlights=options['highlights'])

    from galaxygetopt.outputfiles import OutputFiles
    off = OutputFiles(name='pause_plot', GGO=opts)
    off.CRR(data=output)
