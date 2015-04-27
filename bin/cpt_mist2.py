#!/usr/bin/env python
import argparse
import tempfile
import shutil
import logging
import os
import math
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()
from Bio import SeqIO
import subprocess
IMAGE_BORDER = 50
IMAGE_BORDER_COORD = '%sx%s' % (IMAGE_BORDER, IMAGE_BORDER)

def add_image_border(outfile, i, j, file_data, img_data, zoom):
    half_height = img_data['orig_height']/2
    half_width = img_data['orig_width']/2

    cmd = ['convert', img_data['orig'],
           '-bordercolor', 'purple',
           '-border', '2x2',
           '-bordercolor', 'gray',
           '-border', IMAGE_BORDER_COORD,
           '-rotate', '-90',
           '-pointsize', '40',
           '-font', 'Ubuntu-Mono-Regular',
           '-fill', 'black', '-annotate',
           '+%s+%s' % (half_height, img_data['orig_width'] + IMAGE_BORDER + 45), file_data[j]['header'],
           '-pointsize', '20',
           ]
    for z in range(0, int(.8 * file_data[j]['seqlen']/zoom), int(BestTick(file_data[j]['seqlen'], 5)/zoom)):
        label = label_formatter(z, zoom)
        cmd += [
            '-fill', 'grey22', '-stroke', 'none', '-annotate', '+%s+%s' % (z + 56, IMAGE_BORDER + img_data['orig_width'] + 20), label
        ]

    cmd += ['-rotate', '90',
            '-pointsize', '40',
            '-fill', 'black', '-annotate',
            '+%s+30' % half_width, file_data[i]['header'],
            '-pointsize', '20',
            '-annotate', '+%s+%s' % (2, IMAGE_BORDER + img_data['orig_height'] + 30),
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
            '-fill', 'grey22', '-stroke', 'none', '-annotate', '+%s+%s' % (z + 56, IMAGE_BORDER), label
        ]

    for z in range(0, int(file_data[i]['seqlen']/zoom), int(BestTick(file_data[i]['seqlen'], 5)/zoom)):
        label = label_formatter(z, zoom)
        cmd += [
            '-fill', 'white', '-stroke', 'grey22', '-strokewidth', '2', '-draw', 'line %s,35 %s,%s' % (z + 51, z + 51, IMAGE_BORDER),
        ]

    for z in range(0, int(file_data[j]['seqlen']/zoom), int(BestTick(file_data[j]['seqlen'], 5)/zoom)):
        cmd += [
            '-fill', 'white', '-stroke', 'grey22', '-strokewidth', '2', '-draw', 'line 35,%s %s,%s' % (z + 51, IMAGE_BORDER, z + 51),
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


def run_gepard(seq_a, seq_b, zoom, output, matrix, window):
    print "Running %s vs %s" % (seq_a, seq_b)
    cmd = ['java', '-jar', '/var/lib/gepard.jar',
           '--seq1', seq_a,
           '--seq2', seq_b,
           '--matrix', matrix,
           '--outfile', output,
           '--word', str(window),
           '--zoom', str(zoom),
           '--silent'
           ]
    subprocess.check_call(cmd)


def resize_image(scale, from_file, to_file):
    cmd = ['convert', '-resize', scale, from_file, to_file]
    subprocess.check_call(cmd)


def mist(file, label, zoom, matrix, window=10, filepath=None):
    inputs = zip(file, label)

    tmpdir = tempfile.mkdtemp(prefix="cpt.mist.")
    # If no path provided, make a directory
    if filepath is None:
        outdir = tempfile.mkdtemp(prefix='cpt.mist.')
    else:
        if os.path.isdir(filepath):
            outdir = filepath
        else:
            try:
                os.mkdir(filepath)
                outdir = filepath
            except:
                # TODO: make this less terrible
                raise Exception("--filepath must be a directory")

    ordering = []
    file_data = {}

    for (f, l) in inputs:
        parsed_new_files = extract_info_from_file(f, l, tmpdir)
        print parsed_new_files
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
                       zoom, gepard_png, matrix, window)
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
           '-border', IMAGE_BORDER_COORD,
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

            if ordering.index(i) == 1:
                cumulative_width += img_array[i][j]['height'] + 2

            if ordering.index(j) == 0:
                cumulative_height += img_array[i][j]['width'] + 2

            gepard_annot_png = os.path.join(tmpdir, 'png', 'a_%s-%s.png' % (i, j))
            add_image_border(gepard_annot_png, i, j, file_data,
                             img_array[i][j], int(zoom))
            img_array[i][j]['annotated_original'] = gepard_annot_png

    # The +1 and +2 are as a result of adding a 1 width purple border, so the border is consistent everywhere.
    convert_arguments_top = []
    convert_arguments_left = []
    left_offset = cumulative_width + IMAGE_BORDER

    current_sum_width = IMAGE_BORDER
    current_sum_height = IMAGE_BORDER

    # Top side
    for i in ordering:
        current_image = img_array[i][ordering[0]]
        convert_arguments_top += [
            '-fill', 'black', '-annotate',
            '+%s+40' % current_sum_width, file_data[i]['header']
        ]
        current_sum_width += current_image['width'] + (2 * inter_image_borders)
        print "CSW: %s" % (current_sum_width)

    # Left side
    for i in ordering:
        current_image = img_array[ordering[0]][i]
        convert_arguments_left += [
            '-fill', 'black', '-annotate',
            '+%s+%s' % (current_sum_height, left_offset), '\n' + file_data[i]['header']
        ]

        current_sum_height += current_image['height'] + (2 * inter_image_borders)
        print "CSH: %s" % (current_sum_height)

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
        '-annotate', '+%s+%s' % (2, current_sum_height + 15),
        (
            "\nProduced by the CPT's MIST (Multiple Interrelated "
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
    print html_page

    for (original, new) in files_to_move:
        shutil.copy(original, outdir)
    print outdir


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Identify shine-dalgarno sequences')
    parser.add_argument('file', nargs='+', type=file, help='Input file')
    parser.add_argument('label', nargs='+', help='Label (per input file)')
    parser.add_argument('--zoom', type=int, help='How zoomed in the image is. Be careful with this option. '
                                       'It represents the number of bases to plot in a single pixel.'
                                       ' For large genomes, this can mean very large images, and '
                                       'should be lowered appropriately. For a value of 49, 50 '
                                       'bases would be considered a single pixel in the output image.'
                                       ' For 1Mbp of genomes totally (say 5 x 200 kb phages), this '
                                       'would result in a 20,000 pixel image.', default=50)
    parser.add_argument('--window', type=int, help='Window size', default=10)
    parser.add_argument('--matrix', choices=['ednaorig.mat', 'pam250.mat',
                                             'edna.mat', 'protidentity.mat',
                                             'blosum62.mat'], help='Comparison Matrix', default='edna.mat')
    parser.add_argument('--filepath', help='Directory for images')
    args = parser.parse_args()

    if len(args.file) != len(args.label):
        raise Exception("Must provide the same number of files and labels")
    result = mist(**vars(args))
