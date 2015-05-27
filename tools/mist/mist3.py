#!/usr/bin/env python
import argparse
import tempfile
import shutil
import os
import math
from string import Template
from Bio import SeqIO
import subprocess

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("mist")


class Misty(object):

    IMAGE_BORDER = 50
    IMAGE_BORDER_COORD = '%sx%s' % (IMAGE_BORDER, IMAGE_BORDER)
    CREDITS = (
        "Produced by the CPT's MISTv3 (Multiple Interrelated Sequence doT plotter). "
        "Written by Eric Rasche <esr\@tamu.edu>.\n"
        "Dot plots produced by the Gepard Dot Plotter by Dr. Jan Krumsiek"
    )

    def __init__(self, window=10, zoom=50, matrix='edna', files_path='mist_images'):
        self.tmpdir = tempfile.mkdtemp(prefix="cpt.mist3.", dir='.')
        self.window = str(window)
        self.zoom = zoom
        self.matrix = matrix
        self.records = []
        self.matrix_data = []

        # Image storage
        self.files_path = files_path
        if not os.path.exists(self.files_path):
            os.makedirs(self.files_path)

    def _get_record(self, record_id):
        matched = [x for x in self.records if x['id'] == record_id]
        if len(matched) == 1:
            return matched[0]
        else:
            raise RuntimeError("%s instances in self.records with ID=%s" % (len(matched), record_id))

    def register_all_files(self, file_list):
        for fasta_file in file_list:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                self.register_record(record)

    def register_record(self, record):
        temp = tempfile.NamedTemporaryFile(dir=self.tmpdir, delete=False)
        file_id = temp.name.rsplit('/')[-1]
        self.records.append({
            'file_path': temp.name,
            'id': file_id,
            'header': record.id,
            'description': record.description,
            'length': len(record.seq),
        })
        SeqIO.write([record], temp, 'fasta')

    def set_matrix(self, matrix):
        # More processing?
        self.matrix_data = matrix

    def generate_matrix(self, mtype='complete'):
        matrix = []
        if mtype == 'complete':
            for i in self.records:
                row = []
                for j in self.records:
                    row.append({
                        'i': i['id'],
                        'j': j['id']
                    })
                matrix.append(row)
        elif mtype == '1vn':
            if len(self.records) < 2:
                raise RuntimeError("1vN not available for fewer than two sequences")
            else:
                row = []
                for j in self.records[1:]:
                    row.append({
                        'i': self.records[0]['id'],
                        'j': j['id']
                    })
                matrix.append(row)
        return matrix


    def run_gepard(self, seq_a, seq_b, output):
        log.info("Running %s vs %s", seq_a, seq_b)
        cmd = ['java', '-jar', '/var/lib/gepard.jar',
               '--seq1', seq_a,
               '--seq2', seq_b,
               '--matrix', self.matrix,
               '--outfile', output,
               '--word', self.window,
               '--zoom', str(self.zoom),
               '--silent'
               ]
        log.debug(' '.join(cmd))
        subprocess.check_call(cmd)

    def add_image_border(self, outfile, i, j, file_data, img_data):
        half_height = img_data['orig_height'] / 2
        half_width = img_data['orig_width'] / 2

        j_range = int(1.0 * file_data[j]['seqlen'] / self.zoom)
        i_range = int(1.0 * file_data[i]['seqlen'] / self.zoom)

        j_ticks = int(self.BestTick(file_data[j]['seqlen'], 5)) / self.zoom
        i_ticks = int(self.BestTick(file_data[i]['seqlen'], 5)) / self.zoom
        # Convert definitions
        GREY_FILL = ['-fill', 'grey22']
        WHITE_FILL = ['-fill', 'white']
        NO_STROKE = ['-stroke', 'none']
        GREY_STROKE = ['-stroke', 'grey22', '-strokewidth', '2']
        GFNS = GREY_FILL + NO_STROKE
        WFGS = WHITE_FILL + GREY_STROKE

        cmd = ['convert', img_data['orig'],
               '-bordercolor', 'purple',
               '-border', '2x2',
               '-bordercolor', 'gray',
               '-border', self.IMAGE_BORDER_COORD,
               '-rotate', '-90',
               '-pointsize', '40',
               '-font', 'Ubuntu-Mono-Regular',
               '-fill', 'black', '-annotate',
               '+%s+%s' % (half_height, img_data['orig_width'] + self.IMAGE_BORDER + 45),
               file_data[j]['header'],
               '-pointsize', '20',
               ]

        for z in range(0, int(0.8 * j_range), j_ticks):
            cmd += GFNS
            cmd += [
                '-annotate',
                '+%s+%s' % (z + 56, self.IMAGE_BORDER + img_data['orig_width'] + 20),
                self.label_formatter(z)
            ]

        cmd += ['-rotate', '90',
                '-pointsize', '40',
                '-fill', 'black', '-annotate',
                '+%s+30' % half_width, file_data[i]['header'],
                '-pointsize', '20',
                '-annotate', '+%s+%s' % (2, self.IMAGE_BORDER + img_data['orig_height'] + 30),
                self.CREDITS
                ]

        for z in range(0, int(.8 * i_range), i_ticks):
            cmd += GFNS
            cmd += [
                '-annotate',
                '+%s+%s' % (z + 56, self.IMAGE_BORDER),
                self.label_formatter(z)
            ]

        for z in range(0, i_range, i_ticks):
            cmd += WFGS
            cmd += [
                '-draw',
                'line %s,35 %s,%s' % (z + 51, z + 51, self.IMAGE_BORDER),
            ]

        for z in range(0, j_range, j_ticks):
            cmd += WFGS
            cmd += [
                '-draw',
                'line 35,%s %s,%s' % (z + 51, self.IMAGE_BORDER, z + 51),
            ]

        cmd.append(outfile)
        log.debug(' '.join(cmd))
        subprocess.check_output(cmd)

    def label_formatter(self, bases):
        label = bases * self.zoom
        if label > 1000000:
            label = "%s Mbp" % int(label / 1000000)
        elif label > 1000:
            label = "%s kbp" % int(label / 1000)
        else:
            label = "%s b" % label
        return label

    @classmethod
    def BestTick(cls, largest, mostticks):
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

    @classmethod
    def resize_image(cls, scale, from_file, to_file):
        cmd = ['convert', '-resize', scale, from_file, to_file]
        log.debug(' '.join(cmd))
        subprocess.check_call(cmd)

    def extract_info_from_file(self, infile):
        ret = []
        for seq_record in SeqIO.parse(infile, "fasta"):
            temp = tempfile.NamedTemporaryFile(dir=self.tmpdir, delete=False)
            file_id = temp.name.rsplit('/')[-1]

            ret.append({
                'fasta_path': temp.name,
                'id': file_id,
                'header': seq_record.id,
                'seqlen': len(seq_record.seq),
            })

            temp.write('>%s\n%s' % (seq_record.id, seq_record.seq))

        return ret

    def mist(self, files):
        ordering = []
        file_data = {}

        for input_file in files:
            parsed_new_files = self.extract_info_from_file(input_file)
            log.info('Parsed out: ' + ','.join(parsed_new_files))
            for f in parsed_new_files:
                ordering.append(f['id'])
                file_data[f['id']] = f

        for subdir in ['png', 'thumb', 'prefinal']:
            directory = os.path.join(self.tmpdir, subdir)
            if not os.path.exists(directory):
                os.makedirs(directory)

        total_seqlen = 0
        for i in file_data:
            total_seqlen += file_data[i]['seqlen']

        # Approx how many pixels wide it is
        approx_size = total_seqlen / float(self.zoom)
        # Goal number of pixels
        rescale = 1000.0 / approx_size
        rescale_p = '%s%%' % (100 * rescale)
        img_array = {}
        for i in ordering:
            img_array[i] = {}
            for j in ordering:
                gepard_png = os.path.join(self.tmpdir, 'png', '%s-%s.png' % (i, j))
                thumb_png = os.path.join(self.tmpdir, 'thumb', '%s-%s.png' % (i, j))
                self.run_gepard(file_data[i]['fasta_path'],
                                file_data[j]['fasta_path'],
                                gepard_png)
                self.resize_image(rescale_p, gepard_png, thumb_png)
                img_array[i][j] = {
                    'orig': gepard_png,
                    'thumb': thumb_png,
                }

        cardinality = len(ordering)
        inter_image_borders = 2

        reordered = []
        for i in ordering:
            for j in ordering:
                reordered.append(os.path.join(self.tmpdir, 'thumb', '%s-%s.png' % (j, i)))
        # Make the monatge
        m0 = os.path.join(self.tmpdir, 'prefinal', 'montage_0.png')
        m1 = os.path.join(self.tmpdir, 'prefinal', 'montage_1.png')
        cmd = ['montage'] + reordered + \
            ['-tile', '%sx%s' % (cardinality, cardinality),
             '-geometry', '+0+0',
             '-border', str(inter_image_borders),
             '-bordercolor', 'purple',
             m0]
        subprocess.check_call(cmd)
        # Add borders and labels
        cmd = [
            'convert', m0,
            '-bordercolor', 'purple',
            '-border', '1x1',
            '-bordercolor', 'gray',
            '-border', self.IMAGE_BORDER_COORD,
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

                gepard_annot_png = os.path.join(self.tmpdir, 'png', 'a_%s-%s.png' % (i, j))
                self.add_image_border(gepard_annot_png, i, j, file_data,
                                      img_array[i][j])
                img_array[i][j]['annotated_original'] = gepard_annot_png

        # The +1 and +2 are as a result of adding a 1 width purple border, so the border is consistent everywhere.
        convert_arguments_top = []
        convert_arguments_left = []
        left_offset = cumulative_width + self.IMAGE_BORDER

        current_sum_width = self.IMAGE_BORDER
        current_sum_height = self.IMAGE_BORDER

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

        cmd = [
            'convert', m1,
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
            self.CREDIS,
            os.path.join(self.tmpdir, 'prefinal', 'large.png')
        ]
        log.debug(' '.join(cmd))
        subprocess.check_output(cmd)

        imagemap = ""

        files_to_move = [
            (os.path.join(self.tmpdir, 'prefinal', 'large.png'), 'large')
        ]

        image_template = '<area shape="rect" coords="{x1},{y1},{x2},{y2}" alt="{alt}" href="{href}" />\n'
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
                imagemap += image_template.format({
                    'x1': int(cur_x),
                    'y1': int(cur_y),
                    'x2': int(cur_x + width + 2),
                    'y2': int(cur_y + height + 2),
                    'alt': "%s vs %s" % (file_data[i]['header'], file_data[j]['header']),
                    'href': 'a_%s-%s.png' % (i, j)
                })

                cur_x += width + 4
            cur_y += tmp_height + 4

        for (original, new) in files_to_move:
            produced_files = of.subCRR(data_format="dummy", format_as="Dummy",
                                    filename=new, extension='png', data="dummy")
            destination = produced_files[0]
            shutil.copy(original, destination)

def mist_wrapper(files, zoom=50, matrix='edna', plot_type='complete', files_path='mist_images'):
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
    """

    m = Misty(window=10, zoom=zoom, matrix=matrix, files_path=files_path)

    if plot_type == '2up' and len(files) != 2 and matrix not in ('protidentity', 'blosum62'):
        idx = 0
        # Pull top two sequences.
        for fasta_file in files:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                m.register_record(record)

                # Exit after we've seen two sequences
                idx += 1
                if idx == 2:
                    break
    else:
        m.register_all_files(files)


    if plot_type == 'complete':
        # ALL sequences are used.
        m.set_matrix(m.generate_matrix(mtype='complete'))
    elif plot_type == '2up':
        # Just two sequences.
        if len(files) == 2 and matrix in ('protidentity', 'blosum62'):
            m.set_matrix(m.generate_matrix(mtype='complete'))
        else:
            # Co-op the 1vn to do a single plot, rather than a "complete" plot
            m.set_matrix(m.generate_matrix(mtype='1vn'))
    elif plot_type == '1vn':
        m.set_matrix(m.generate_matrix(mtype='1vn'))
    else:
        raise ValueError("Unknown plot type %s" % plot_type)

    # image_map will be returned from MIST
    image_map = ""

    return html_page % image_map


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MIST v3', epilog='')
    parser.add_argument('files', nargs='+', type=file, help='Fasta sequences')

    parser.add_argument('--zoom', type=int, help='# bases / pixel', default=50)
    parser.add_argument('--matrix', type=str, choices=['ednaorig', 'pam250',
                                                       'edna', 'protidentity',
                                                       'blosum62'],
                        help='# bases / pixel', default='edna')

    parser.add_argument('--plot_type', choices=['2up', '1vn', 'complete'], help='Plot type', default='complete')

    parser.add_argument('--files_path', type=str, help='Image directory', default='mist_images')

    args = parser.parse_args()
    print mist_wrapper(**vars(args))
