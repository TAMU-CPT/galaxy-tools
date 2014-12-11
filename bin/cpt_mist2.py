#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
from galaxygetopt.outputfiles import OutputFiles
import tempfile
import shutil
import logging
import os
import math
logging.basicConfig(level=logging.INFO)
from Bio import SeqIO
import subprocess


__doc__ = """
Multiple Interrlated Sequence doTplotter
========================================

"""


def add_image_border(outfile, i, j, file_data, img_data, zoom):
    half_height = img_data['orig_height']/2
    half_width = img_data['orig_width']/2

    cmd = ['convert', img_data['orig'],
           '-bordercolor', 'purple',
           '-border', '2x2',
           '-bordercolor', 'gray',
           '-border', '50x50',
           '-rotate', '-90',
           '-pointsize', '40',
           '-font', 'Ubuntu-Mono-Regular',
           '-fill', 'black', '-annotate',
           '+%s+%s' % (half_height, img_data['orig_width'] + 50 + 45), file_data[j]['header'],
           '-pointsize', '20',
           ]
    for z in range(0, int(.8 * file_data[j]['seqlen']/zoom), int(BestTick(file_data[j]['seqlen'], 5)/zoom)):
        label = label_formatter(z, zoom)
        cmd += [
            '-fill', 'grey22', '-stroke', 'none', '-annotate', '+%s+%s' % (z + 56, 50 + img_data['orig_width'] + 20), label
        ]

    cmd += ['-rotate', '90',
            '-pointsize', '40',
            '-fill', 'black', '-annotate',
            '+%s+30' % half_width, file_data[i]['header'],
            '-pointsize', '20',
            '-annotate', '+%s+%s' % (2, 50 + img_data['orig_height'] + 30),
            (
                "Produced by the CPT's MIST (Multiple Interrelated "
                "Sequence doT plotter). "
                "Written by Eric Rasche <rasche.eric\@yandex.ru>.\n"
                "Dot plots produced by the Gepard Dot Plotter by Dr. Jan Krumsiek"
            ),
            ]
    for z in range(0, int(.8 * file_data[i]['seqlen']/zoom), int(BestTick(file_data[i]['seqlen'], 5)/zoom)):
        label = label_formatter(z, zoom)
        cmd += [
            '-fill', 'grey22', '-stroke', 'none', '-annotate', '+%s+50' % (z + 56), label
        ]

    for z in range(0, int(file_data[i]['seqlen']/zoom), int(BestTick(file_data[i]['seqlen'], 5)/zoom)):
        label = label_formatter(z, zoom)
        cmd += [
            '-fill', 'white', '-stroke', 'grey22', '-strokewidth', '2', '-draw', 'line %s,35 %s,50' % (z + 51, z + 51),
        ]

    for z in range(0, int(file_data[j]['seqlen']/zoom), int(BestTick(file_data[j]['seqlen'], 5)/zoom)):
        cmd += [
            '-fill', 'white', '-stroke', 'grey22', '-strokewidth', '2', '-draw', 'line 35,%s 50,%s' % (z + 51, z + 51),
        ]

    cmd.append(outfile)
    subprocess.check_output(cmd)


def label_formatter(bases, zoom):
    label = bases * zoom
    if label > 1000000:
        label = "%s Mbp" % int(label/1000000)
    elif label > 1000:
        label = "%s kbp" % int(label/1000)
    else:
        label = "%s b" % label
    return label


def BestTick(largest, mostticks):
    # http://stackoverflow.com/a/361687
    minimum = largest / mostticks
    magnitude = 10 ** math.floor(math.log(minimum) / math.log(10))
    residual = minimum / magnitude
    if residual > 5:
        tick = 10 * magnitude
    elif residual > 2:
        tick = 5 * magnitude
    elif residual > 1:
        tick = 2 * magnitude
    else:
        tick = magnitude
    return tick


def extract_info_from_file(infile, label, tmpdir):
    ret = []
    for seq_record in SeqIO.parse(infile, "fasta"):
        temp = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
        file_id = temp.name.rsplit('/')[-1]

        ret.append({
            'fasta_path': temp.name,
            'id': file_id,
            'header': seq_record.id,
            'seqlen': len(seq_record.seq),
        })

        temp.write('>%s\n%s' % (seq_record.id, seq_record.seq))
    # We overwrite the first one with label, if there are multiple, those
    # retain their own labels.
    if len(ret) > 0:
        ret[0]['header'] = label
    return ret


def run_gepard(seq_a, seq_b, zoom, output, matrix):
    print "Running %s vs %s" % (seq_a, seq_b)
    cmd = ['java', '-jar', '/var/lib/gepard.jar',
           '--seq1', seq_a,
           '--seq2', seq_b,
           '--matrix', matrix,
           '--outfile', output,
           '--zoom', zoom,
           '--silent'
           ]
    subprocess.check_call(cmd)


def resize_image(scale, from_file, to_file):
    cmd = ['convert', '-resize', scale, from_file, to_file]
    subprocess.check_call(cmd)


def mist(ggo, file, label, zoom, matrix, *args, **kwargs):
    inputs = zip(file, label)

    tmpdir = tempfile.mkdtemp(prefix="cpt.mist.")

    ordering = []
    file_data = {}

    for input_file in inputs:
        parsed_new_files = extract_info_from_file(input_file[0], input_file[1],
                                                  tmpdir)
        for f in parsed_new_files:
            ordering.append(f['id'])
            file_data[f['id']] = f

    for subdir in ['png', 'thumb', 'prefinal']:
        directory = os.path.join(tmpdir, subdir)
        if not os.path.exists(directory):
            os.makedirs(directory)

    total_seqlen = 0
    for i in file_data:
        total_seqlen += file_data[i]['seqlen']

    # Approx how many pixels wide it is
    approx_size = total_seqlen/float(zoom)
    # Goal number of pixels
    rescale = 1000.0 / approx_size
    rescale_p = '%s%%' % (100 * rescale)
    img_array = {}
    for i in ordering:
        img_array[i] = {}
        for j in ordering:
            gepard_png = os.path.join(tmpdir, 'png', '%s-%s.png' % (i, j))
            thumb_png = os.path.join(tmpdir, 'thumb', '%s-%s.png' % (i, j))
            run_gepard(file_data[i]['fasta_path'], file_data[j]['fasta_path'],
                       zoom, gepard_png, matrix)
            resize_image(rescale_p, gepard_png, thumb_png)
            img_array[i][j] = {
                'orig': gepard_png,
                'thumb': thumb_png,
            }

    cardinality = len(ordering)
    inter_image_borders = 2

    reordered = []
    for i in ordering:
        for j in ordering:
            reordered.append(os.path.join(tmpdir,
                                          'thumb', '%s-%s.png' % (j, i)))
    # Make the monatge
    m0 = os.path.join(tmpdir, 'prefinal', 'montage_0.png')
    m1 = os.path.join(tmpdir, 'prefinal', 'montage_1.png')
    cmd = ['montage'] + reordered + \
        ['-tile', '%sx%s' % (cardinality, cardinality),
         '-geometry', '+0+0',
         '-border', str(inter_image_borders),
         '-bordercolor', 'purple',
         m0]
    subprocess.check_call(cmd)
    # Add borders and labels
    cmd = ['convert', m0,
           '-bordercolor', 'purple',
           '-border', '1x1',
           '-bordercolor', 'gray',
           '-border', '50x50',
           m1
           ]
    subprocess.check_call(cmd)

    cumulative_width = 0
    cumulative_height = 0
    for i in ordering:
        for j in ordering:
            current_image = img_array[i][j]
            output = subprocess.check_output(['identify',
                                              current_image['orig']])
            size = output.split(' ')[3]
            (w, h) = size[0:size.index('+')].split('x')
            img_array[i][j]['width'] = int(w) * rescale
            img_array[i][j]['height'] = int(h) * rescale
            img_array[i][j]['orig_width'] = int(w)
            img_array[i][j]['orig_height'] = int(h)

            if ordering.index(i) == 0:
                cumulative_width += img_array[i][j]['height'] + 2

            if ordering.index(j) == 0:
                cumulative_height += img_array[i][j]['width'] + 2

            gepard_annot_png = os.path.join(tmpdir, 'png', 'a_%s-%s.png' % (i, j))
            add_image_border(gepard_annot_png, i, j, file_data,
                             img_array[i][j], int(zoom))
            img_array[i][j]['annotated_original'] = gepard_annot_png

    # The +1 and +2 are as a result of adding a 1 width purple border, so the border is consistent everywhere.
    current_sum = 51 + inter_image_borders
    convert_arguments_top = []
    convert_arguments_left = []
    left_offset = cumulative_width + 51 + 30

    for i in ordering:
        current_image = img_array[i][i]
        convert_arguments_top += [
            '-fill', 'black', '-annotate',
            '+%s+40' % current_sum, file_data[i]['header']
        ]
        print "CS: %s LO: %s" % (current_sum, left_offset)
        convert_arguments_left += [
            '-fill', 'black', '-annotate',
            '+%s+%s' % (current_sum - 10, left_offset), file_data[i]['header']
        ]

        current_sum += current_image['width'] + (2 * inter_image_borders)

    cmd = ['convert', m1,
           '-rotate', '-90',
           '-pointsize', '20',
           '-font', 'Ubuntu-Mono-Regular',
           ]
    cmd += convert_arguments_left
    cmd += ['-rotate', '90']
    cmd += convert_arguments_top
    cmd += [
        '-pointsize', '14',
        '-annotate', '+%s+%s' % (2, left_offset),
        (
            "Produced by the CPT's MIST (Multiple Interrelated "
            "Sequence doT plotter). "
            "Written by Eric Rasche <rasche.eric\@yandex.ru>.\n"
            "Dot plots produced by the Gepard Dot Plotter by Dr. Jan Krumsiek"
        ),
        os.path.join(tmpdir, 'prefinal', 'large.png')
    ]
    print ' '.join(cmd)
    subprocess.check_output(cmd)

    imagemap = ""

    files_to_move = [
        (os.path.join(tmpdir, 'prefinal', 'large.png'), 'large')
    ]

    cur_y = 51
    for j in ordering:
        cur_x = 51
        tmp_height = 0
        for i in ordering:
            current = img_array[i][j]
            width = current['width']
            height = current['height']
            tmp_height = height
            files_to_move.append((img_array[i][j]['annotated_original'],
                                  'a_%s-%s' % (i, j)))
            print current
            imagemap += '<area shape="rect" coords="%s,%s,%s,%s" alt="%s" href="%s" />\n' \
                % (int(cur_x),
                   int(cur_y),
                   int(cur_x+width+2),
                   int(cur_y+height+2),
                    "%s vs %s" % (file_data[i]['header'],
                                  file_data[j]['header']),
                   'a_%s-%s.png' % (i, j)
                   )
            cur_x += width + 4
        cur_y += tmp_height + 4


    html_page = """
    <!DOCTYPE html>
    <html>
    <body>
    <h1>Mist Results</h1>
    <p>Each section of mist output is now clickable to view a higher resolution version of that subsection</p>
    <img src="large.png" usemap="#mistmap">
    <map name="mistmap">
    %s
    </map>
    </body>
    </html>
    """ % imagemap

    of = OutputFiles(name='dotplot', GGO=ggo)
    of.CRR(data=html_page)

    for (original, new) in files_to_move:
        produced_files = of.subCRR(data_format="dummy", format_as="Dummy",
                                   filename=new, extension='png', data="dummy")
        destination = produced_files[0]
        shutil.copy(original, destination)


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    opts = GGO(
        options=[
            ['file', 'Input file',
             {'required': True, 'multiple': True, 'validate': 'File/Input'}],
            ['label', 'Label (per input file)',
             {'required': True, 'multiple': True, 'validate': 'String'}],
            ['zoom', 'How zoomed in the image is. Be careful with this option. It represents the number of bases to plot in a single pixel. For large genomes, this can mean very large images, and should be lowered appropriately. For a value of 50, 50 bases would be considered a single pixel in the output image. For 1Mbp of genomes totaly (say 5 x 200 kb phages), this would result in a 20,000 pixel image.',
             {'required': True, 'validate': 'String', 'default': 50}],
            ['matrix', 'Comparison Matrix',
             {'required': True, 'validate': 'Option', 'options': {
                 'ednaorig.mat': 'Extended DNA (Original)',
                 'pam250.mat': 'Pam 250',
                 'edna.mat': 'Extended DNA',
                 'protidentity.mat': 'Protein Identity',
                 'blosum62': 'Blosum62',
             }, 'default': 'edna.mat'}]
        ],
        outputs=[
            [
                'dotplot',
                'MIST Plot',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'extracted',
                    'data_format': 'text/html',
                    'default_format': 'HTML',
                }
            ],
        ],
        defaults={
            'appid': 'edu.tamu.cpt.mist',
            'appname': 'MIST',
            'appvers': '0.1',
            'appdesc': 'Multiple Interrelated Sequence doT plotter. Uses a stripped down vesion of Gepard (Dr. Jan Krumsiek/HelmholtzZentrum/IBIS) for dot plotting.',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()
    result = mist(ggo=opts, **options)
    #svg_xml = generate_viz(blast_file=options['blast'])

