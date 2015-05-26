#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import logging
logging.basicConfig(level=logging.INFO)


__doc__ = """
Blast Alignment Viz.
====================

"""


def coallesce(blast_hits=[], query=None, subject=None):
    # In a set of blast hits, we order based on left to right in our query
    # genome
    # We need to collapse these into blocks (generally 1 or 2, depending on
    # re-opening)
    blocks = []
    current_block = None
    overlap_fraction = int(subject * 0.01)
    # There's no basis for this number
    gap_fraction = 2 * overlap_fraction
    # map of subject start/ends to query sequence locations
    start_end_map = {}
    print "Q:%s S:%s" % (query, subject)
    print "Allowing %s NTs overlap" % overlap_fraction
    for location in sorted(blast_hits, key=lambda row: int(row[0])):
        start_end_map[location[2]] = location[0]
        start_end_map[location[3]] = location[1]
        #print location
        if current_block is None:
            current_block = location[2:4]
            # Increasing
        else:
            if location[2] > current_block[1] - overlap_fraction and \
               location[2] < current_block[1] + gap_fraction:
                #print "Increasing"
                # Extend
                current_block[1] = location[2]
            elif location[2] > current_block[1] - overlap_fraction and \
                    location[2] >= current_block[1] + gap_fraction:
                #print "Gap"
                blocks.append(current_block)
                current_block = None
            else:
                #print "Decreasing"
                blocks.append(current_block)
                current_block = None
        current_block = location[2:4]
    blocks.append(current_block)

    location_map = []
    for block in blocks:
        location_map.append([
            (start_end_map[block[0]], start_end_map[block[1]]),
            (block[0], block[1]),
        ])
    return location_map


def generate_viz(blast_file=None):
    # Unique list of query titles, for use in output
    queries = {}
    length_database = {}
    for line in blast_file.readlines():
        data = line.split('\t')
        if data[0] not in queries:
            queries[data[0]] = {}
        if data[1] not in queries[data[0]]:
            queries[data[0]][data[1]] = []
        queries[data[0]][data[1]].append([float(x) for x in [
            data[6], data[7], data[8], data[9], data[10], data[11]
        ]])
        length_database[data[0]] = int(data[22])
        length_database[data[1]] = int(data[23])

    import svgwrite
    svg_document = svgwrite.Drawing(size=("1200px", "1200px"))

    current_x = 20
    current_y = 20

    for genome in queries:
        svg_document.add(svg_document.text(genome, insert=(current_x, current_y)))
        current_y += 30
        for hit in queries[genome]:
            svg_document.add(svg_document.text(queries, insert=(current_x+20, current_y)))
            current_y += 30
            hits = coallesce(blast_hits=queries[genome][hit],
                             query=length_database[genome],
                             subject=length_database[hit])
            svg_document.add(svg_document.rect(
                insert=(100, current_y),
                size=("%spx" % 1000, "%spx" % 20),
                stroke_width="1",
                stroke="black",
                fill="red"))
            current_y += 20

            for box in hits:
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
            ['blast', 'Blast file to plot',
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
            'appid': 'edu.tamu.cpt.genbank.blastalignment',
            'appname': 'Blast Alignment Visualization',
            'appvers': '1.94',
            'appdesc': 'plots blast alignments',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    svg_xml = generate_viz(blast_file=options['blast'])

    from galaxygetopt.outputfiles import OutputFiles
    of = OutputFiles(name='plot', GGO=opts)
    of.CRR(data=svg_xml)
