#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import tempfile
import logging
import os
logging.basicConfig(level=logging.INFO)
from Bio import SeqIO
import subprocess


__doc__ = """
Multiple Interrlated Sequence doTplotter
========================================

"""


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


def run_gepard(seq_a, seq_b, zoom, output):
    cmd = ['java', '-jar', 'gepard.jar',
           '--seq1', seq_a,
           '--seq2', seq_b,
           '--matrix', 'edna.mat',
           '--outfile', output,
           '--zoom', zoom
           ]
    print ' '.join(cmd)
    subprocess.check_call(cmd)


def resize_image(scale, from_file, to_file):
    cmd = ['convert', '-resize', scale, from_file, to_file]
    subprocess.check_call(cmd)


def mist(file, label, zoom, java_7_bin, *args, **kwargs):
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

    for subdir in ['png', 'thumb', 'prefinal', 'final']:
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
                       zoom, gepard_png)
            resize_image(rescale_p, gepard_png, thumb_png)
            img_array[i][j] = {
                'from': i,
                'to': j,
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
                                              current_image['thumb']])
            size = output.split(' ')[3]
            (w, h) = size[0:size.index('+')].split('x')
            img_array[i][j]['width'] = int(w)
            img_array[i][j]['height'] = int(h)

            if ordering.index(i) == 0:
                cumulative_width += int(w)

            if ordering.index(j) == 0:
                cumulative_height += int(h)

    #import pprint; pprint.pprint(img_array)
    border_size = 2 * inter_image_borders * cardinality
    # The +1 and +2 are as a result of adding a 1 width purple border, so the border is consistent everywhere.
    total_size_one_border = cumulative_width + border_size
    current_sum = 1 + inter_image_borders
    convert_arguments_top = []
    convert_arguments_left = []
    left_offset = total_size_one_border + 20

    # Original, bad
    # -pointsize 24 -font Ubuntu-Mono-Regular -fill black -annotate +3+1385 ANI
    # -fill black -annotate +458+1385 Karma -fill black -annotate +913+1385 fasta
    # -rotate 90 -fill black -annotate +3+30 ANI -fill black -annotate +458+30 Karma
    # -fill black -annotate +913+30 fasta -pointsize 14 -annotate +2+1385 "test"

    # Manual, good
    # -pointsize 20 -font Ubuntu-Mono-Regular -fill black -annotate +53+1085 ANI
    # -fill black -annotate +508+1085 Karma -fill black -annotate +963+1085 fasta
    # -rotate 90 -fill black -annotate +53+40 ANI -fill black -annotate +508+40 Karma
    # -fill black -annotate +963+40 fasta -pointsize 14 -annotate +2+1100 "test"

    for i in ordering:
        current_image = img_array[i][i]
        print current_image
        convert_arguments_top += [
            '-fill', 'black', '-annotate',
            '+%s+30' % current_sum, file_data[i]['header']
        ]
        convert_arguments_left += [
            '-fill', 'black', '-annotate',
            '+%s+%s' % (current_sum, left_offset), file_data[i]['header']
        ]

        current_sum += current_image['width'] + (2 * inter_image_borders)

    m2 = os.path.join(tmpdir, 'prefinal', 'montage_2.png')
    cmd = ['convert', m1,
           '-rotate', '-90',
           '-pointsize', '24',
           '-font', 'Ubuntu-Mono-Regular',
           ]
    cmd += convert_arguments_left
    cmd += ['-rotate', '90']
    cmd += convert_arguments_top
    cmd += [
        '-pointsize', '14',
        '-annotate', '+%s+%s' % (2, left_offset),
        "Produced by the CPT's MIST (Multiple Interrelated Sequence doT plotter). Written by Eric Rasche <rasche.eric\@yandex.ru>.\nDot plots produced by the Gepard Dot Plotter by Dr. Jan Krumsiek",
        os.path.join(tmpdir, 'final', 'large.png')
    ]
    print ' '.join(cmd)
    subprocess.check_output(cmd)


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
            ['java_7_bin', 'Location of the Java binary',
             {'required': True, 'default': '/usr/lib/jvm/java-7-openjdk-amd64/bin/java', 'hidden': True, 'validate': 'String'}],
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
    result = mist(**options)
    #svg_xml = generate_viz(blast_file=options['blast'])

    #from galaxygetopt.outputfiles import OutputFiles
    #of = OutputFiles(name='plot', GGO=opts)
    #of.CRR(data=svg_xml)
