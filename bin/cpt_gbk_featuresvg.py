#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import logging
logging.basicConfig(level=logging.INFO)

__doc__ = """
Print Features as Colored Boxes
===============================

c.f. Alan Davidson's software
"""


def get_boxes(features):
    boxes = []
    for feature in features:
        if 'color' in feature.qualifiers:
            boxes.append(feature.qualifiers['color'])
        elif 'colour' in feature.qualifiers:
            boxes.append(feature.qualifiers['colour'])
        else:
            boxes.append(None)
    return boxes


def plot_boxes(parent=None, box_size=20):
    from Bio import SeqIO
    records = list(SeqIO.parse(parent, "genbank"))
    extracted = []
    name_list = []
    max_len = 0
    max_size = len(records)
    for i in range(len(records)):
        # Only do CDSs
        next_row = get_boxes([f for f in records[i].features if f.type ==
                              'CDS'])
        extracted.append(next_row)
        if len(next_row) > max_len:
            max_len = len(next_row)
        desc = records[i].description
        if ',' in desc:
            desc = desc[:desc.index(',')]
        name_list.append(records[i].id + " " + desc)
    import svgwrite

    left_side_gap = 400

    # Convenience functions to calculate positioning
    def calc_x(box_idx):
        return (box_size+1) * box_idx + left_side_gap

    def calc_y(box_idx):
        return (box_size+5) * box_idx

    # Box_size + 1 accounts for border on one side, then +1 for last border
    width = calc_x(max_len) + 1
    # Height is calculated much the same, plus some spacing (4px)
    height = calc_y(max_size) + 1
    svg_document = svgwrite.Drawing(size=("%spx" % width, "%spx" % height))

    for row_idx, gbk_name in zip(range(len(extracted)), name_list):
        row = extracted[row_idx]
        svg_document.add(svg_document.text(gbk_name, insert=(calc_x(0) - left_side_gap, calc_y(row_idx+0.5))))
        for box_idx in range(len(row)):
            if row[box_idx] is not None:
                # [0] is so we always take the "First" specified colour
                rgb = row[box_idx][0].split(' ')
                if len(rgb) == 3:
                    color = "rgb(%s)" % ','.join(rgb)
                else:
                    print rgb
            else:
                if row_idx % 2 == 0:
                    color = "rgb(220,220,220)"
                else:
                    color = "rgb(250,250,250)"
            svg_document.add(svg_document.rect(
                insert=(calc_x(box_idx), calc_y(row_idx)),
                size=("%spx" % box_size, "%spx" % box_size),
                stroke_width="1",
                stroke="black",
                fill=color))

    return svg_document.tostring()


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['gbk', 'Genbank file to plot',
             {'required': True, 'validate': 'File/Input'}],
        ],
        outputs=[
            [
                'plot',
                'SVG Plot',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'extracted',
                    'data_format': 'text/plain',
                    'default_format': 'TXT',
                }
            ],
        ],
        defaults={
            'appid': 'edu.tamu.cpt.genbank.DavidsonPlot',
            'appname': 'Genbank Feature Array Plot',
            'appvers': '1.94',
            'appdesc': 'plots features as an array of colored boxes',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    svg_xml = plot_boxes(parent=options['gbk'])

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='plot', GGO=opts)
    of.CRR(data=svg_xml)
