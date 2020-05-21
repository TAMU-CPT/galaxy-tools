#!/usr/bin/env python 
import argparse
import time
import tempfile
import pprint
import shutil
import os
import math
from string import Template
from Bio import SeqIO
import subprocess

import logging

FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s:%(funcName)s()] %(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT)
log = logging.getLogger("mist")
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

MONTAGE_BORDER = 50
IMAGE_BORDER = 2

MONTAGE_BORDER_COORD = "%sx%s" % (MONTAGE_BORDER, MONTAGE_BORDER)
IMAGE_BORDER_COORD = "%sx%s" % (IMAGE_BORDER, IMAGE_BORDER)

DOUBLE_IMAGE_BORDER_COORD = "%sx%s" % (2 * IMAGE_BORDER, 2 * IMAGE_BORDER)

MONTAGE_BORDER_COLOUR = "grey70"
IMAGE_BORDER_COLOUR = "purple"
LABEL_COLOUR = "grey22"
TICK_LENGTH = 0.2 * MONTAGE_BORDER

TYPEFONT = "Ubuntu-Mono"

CREDITS = (
    "CPT's MISTv3\n"
    "GPLv3 (C) 2015 Helena Rasche <esr\@tamu.edu>\n"
    "Dot plots by Gepard Dot Plotter by Dr. Jan Krumsiek"
)


class FancyRecord(object):
    def __init__(self, record, tmpdir):
        self.temp = tempfile.NamedTemporaryFile(mode='w', dir=tmpdir, delete=False, suffix=".fa")
        self.temp_path = self.temp.name
        self.id = self.temp_path.rsplit("/")[-1]
        self.record = record
        self.header = record.id

        # Description include record.id
        sd = record.description.strip().split(" ")
        if len(sd) > 1:
            self.description = " ".join(sd[1:])
        else:
            self.description = ""

        self.length = len(record.seq)

        # Serialize to disk
        self._write(self.temp)
        self.temp.flush()
        self.temp.close()

    def _write(self, handle):
        SeqIO.write([self.record], handle, "fasta")


class Subplot(object):
    def __init__(self, i, j, tmpdir, zoom):
        self.i = i
        self.j = j
        # Mist information
        self.tmpdir = tmpdir
        self.zoom = zoom

        self.original_path = None
        self.thumb_path = None
        self.annotated_original_path = None

    def safe(self, string):
        return "".join(
            [
                c
                for c in string.strip()
                if c.isalpha() or c.isdigit() or c == " " or c in ("[", "]", "-", "_")
            ]
        ).rstrip()

    def move_to_dir(self, files_path):
        """Move a specific image that's part of a subplot to an output
        directory and return the name
        """
        destination_fn = (
            self.safe("%s_vs_%s_[annotated]" % (self.i.header, self.j.header)) + ".png"
        )
        destination = os.path.join(files_path, destination_fn)
        log.debug("cp %s %s", self.annotated_original_path, destination)
        if (
            self.annotated_original_path is not None
            and os.path.exists(self.annotated_original_path)
            and not os.path.exists(destination)
        ):
            shutil.move(self.annotated_original_path, destination)
        return destination_fn

    def get_description(self):
        return "%s v %s" % (self.i.header, self.j.header)

    def __repr__(self):
        return "Subplot [%s]" % self.get_description()

    def run_gepard(self, matrix, window, global_rescale="35%"):
        """Run gepard on two sequences, with a specified output file
        """
        log.info("Running Gepard on %s", self.get_description())

        destination_fn = (
            self.safe("%s_vs_%s_orig" % (self.i.header, self.j.header)) + ".png"
        )
        self.original_path = os.path.join(self.tmpdir, destination_fn)

        cmd = [
            "java",
            "-jar",
            os.path.join(SCRIPT_DIR, "gepard.jar"),
            "--seq1",
            self.j.temp_path,
            "--seq2",
            self.i.temp_path,
            "--matrix",
            matrix + ".mat",
            "--outfile",
            self.original_path,
            "--word",
            str(window),
            "--zoom",
            str(self.zoom),
            "--silent",
        ]
        log.debug(subprocess.list2cmdline(cmd))

        failure_count = 0
        while True:
            try:
                subprocess.check_call(cmd)
                break
            except subprocess.CalledProcessError:
                failure_count += 1
                log.warn("sp.CPE FC %s", failure_count)
                if failure_count > 3:
                    break
                time.sleep(1)

        # Generate border/individualised images
        log.debug("    Annotating")
        destination_fn = (
            self.safe("%s_vs_%s_annotated" % (self.i.header, self.j.header)) + ".png"
        )
        self.annotated_original_path = os.path.join(self.tmpdir, destination_fn)
        self.annotate_image(self.original_path, self.annotated_original_path)

        # Generate thumbnail version of base image
        log.debug("    Resizing")
        destination_fn = (
            self.safe("%s_vs_%s_thumb" % (self.i.header, self.j.header)) + ".png"
        )
        self.thumb_path = os.path.join(self.tmpdir, destination_fn)
        Misty.resize_image(global_rescale, self.original_path, self.thumb_path)

    def get_thumb_dims(self):
        """
        For NxN sized images, this avoids running `identify` when it won't be
        used (e.g i=3,j=4)
        """
        if not hasattr(self, "thumb_dims") or self.thumb_dims is None:
            self.thumb_dims = Misty.obtain_image_dimensions(self.thumb_path)
        return self.thumb_dims

    def label_formatter(self, bases):
        if bases > 1000000:
            label = "%s Mbp" % int(bases / 1000000)
        elif bases > 1000:
            label = "%s kbp" % int(bases / 1000)
        else:
            label = "%s b" % bases
        return label

    def annotate_image(self, infile, outfile):
        original_dims = Misty.obtain_image_dimensions(infile)

        half_height = (original_dims[1] / 2) + MONTAGE_BORDER + IMAGE_BORDER
        half_width = (original_dims[0] / 2) + MONTAGE_BORDER + IMAGE_BORDER

        def char_width(font_size):
            """approximate pixel width of a single character at Xpt font

            Assumes that 40pt font results in 20px characters
            """
            return int(float(font_size) / 2)

        def char_height(font_size):
            """approximate pixel height of a single character at Xpt font

            Assumes that 40pt font results in 30px characters
            """
            return int(float(font_size) * 30 / 40)

        def est_pixels(string, font_size):
            """guess pixel width of a string at a given font size
            """
            return char_width(font_size) * len(string)

        j_ticks = int(Misty.BestTick(self.j.length, 5))
        i_ticks = int(Misty.BestTick(self.i.length, 5))

        # Convert definitions
        GREY_FILL = ["-fill", LABEL_COLOUR]
        NO_FILL = ["-fill", "none"]
        NO_STROKE = ["-stroke", "none"]
        GREY_STROKE = ["-stroke", LABEL_COLOUR, "-strokewidth", "2"]
        GFNS = GREY_FILL + NO_STROKE
        NFGS = NO_FILL + GREY_STROKE
        # Font for labels
        FONT_SPEC = GFNS + ["-font", TYPEFONT]
        FONT_10pt = FONT_SPEC + ["-pointsize", "10"]
        FONT_20pt = FONT_SPEC + ["-pointsize", "20"]
        FONT_30pt = FONT_SPEC + ["-pointsize", "30"]

        cmd = [
            "convert",
            infile,
            "-bordercolor",
            IMAGE_BORDER_COLOUR,
            "-border",
            IMAGE_BORDER_COORD,
            "-bordercolor",
            MONTAGE_BORDER_COLOUR,
            "-border",
            MONTAGE_BORDER_COORD,
        ]

        # Rotate -90, apply row header at bottom.
        primary_header = self.i.header
        secondary_head = self.i.description
        cmd += (
            ["-rotate", "-90",]
            + FONT_30pt
            + [
                # Side label (i/row)
                "-annotate",
                "+%s+%s"
                % (
                    half_height - 0.5 * est_pixels(self.i.header, 30),
                    original_dims[0]
                    + MONTAGE_BORDER
                    + 2 * IMAGE_BORDER
                    + char_height(30 + 20)
                    # 30 = primary header
                    # 20 = tick labels
                ),
                primary_header,
            ]
        )

        if est_pixels(self.i.description, 10) < original_dims[1]:
            cmd += FONT_10pt + [
                # Side label (i/row)
                "-annotate",
                "+%s+%s"
                % (
                    half_height - 0.5 * est_pixels(self.i.description, 10),
                    original_dims[0]
                    + MONTAGE_BORDER
                    + 2 * IMAGE_BORDER
                    + char_height(30 + 20 + 10 + 4)
                    # 30 = primary header
                    # 20 = tick labels
                    # 10 = secondary header height
                    # 4  = line spacing
                ),
                secondary_head,
            ]

        # Apply row ticks labels at bottom
        for z in range(0, self.i.length, i_ticks):

            image_side_percentage = float(z) / self.i.length
            x = (
                (image_side_percentage * original_dims[1])
                + MONTAGE_BORDER
                + IMAGE_BORDER
            )
            y = MONTAGE_BORDER + original_dims[0] + (2 * IMAGE_BORDER)

            # Apply ticks
            cmd += NFGS
            cmd += [
                "-draw",
                "line %s,%s %s,%s" % (x, y, x, y + TICK_LENGTH),
            ]

            # Keep text from running off the edge.
            space_to_end_of_image = (1 - image_side_percentage) * original_dims[1]
            if space_to_end_of_image - est_pixels(self.label_formatter(z), 20) < 0:
                continue

            # Text label
            cmd += FONT_20pt
            cmd += ["-annotate", "+%s+%s" % (x + 5, y + 15), self.label_formatter(z)]

        # Rotate back to final rotation
        primary_header = self.j.header
        secondary_head = self.j.description
        cmd += (
            [
                "-rotate",
                "90",
                # Top label (j/column)
            ]
            + FONT_30pt
            + [
                "-annotate",
                "+%s+%s"
                % (
                    half_width - 0.5 * est_pixels(self.j.header, 30),
                    MONTAGE_BORDER - char_height(20 + 10 + 8),
                ),
                primary_header,
            ]
            + FONT_10pt
            + [
                # Credits
                "-annotate",
                "+%s+%s" % (2, MONTAGE_BORDER + original_dims[1] + 2 * IMAGE_BORDER),
                "\n" + CREDITS,
            ]
        )

        if est_pixels(self.j.description, 10) < original_dims[0]:
            cmd += FONT_10pt + [
                "-annotate",
                "+%s+%s"
                % (
                    half_width - 0.5 * est_pixels(self.j.description, 10),
                    MONTAGE_BORDER - char_height(20 + 4)
                    # 4  = line spacing
                ),
                secondary_head,
            ]

        # Apply col ticks along top
        for z in range(0, self.j.length, j_ticks):
            image_side_percentage = float(z) / self.j.length
            x = (
                (image_side_percentage * original_dims[0])
                + MONTAGE_BORDER
                + IMAGE_BORDER
            )
            y = MONTAGE_BORDER - 1

            # Ticks
            cmd += NFGS
            cmd += [
                "-draw",
                "line %s,%s %s,%s" % (x, y, x, y - TICK_LENGTH),
            ]

            # Keep text from running off the edge.
            space_to_end_of_image = (1 - image_side_percentage) * original_dims[0]
            if space_to_end_of_image - est_pixels(self.label_formatter(z), 20) < 0:
                continue

            # Text labels
            cmd += FONT_20pt
            cmd += ["-annotate", "+%s+%s" % (x + 5, y), self.label_formatter(z)]

        cmd.append(outfile)
        log.debug(subprocess.list2cmdline(cmd))
        subprocess.check_output(cmd)


class Misty(object):
    """MIST Class for building MIST Plots
    """

    def __init__(self, window=10, zoom=50, matrix="edna", files_path="mist_images"):
        self.tmpdir = tempfile.mkdtemp(prefix="cpt.mist3.", dir=".")
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
        for i, record in enumerate(self.records):
            if record.id == record_id:
                return record

    def _get_record_idx(self, record_id):
        for i, record in enumerate(self.records):
            if record.id == record_id:
                return i

        raise RuntimeError("Could not find record ID=%s" % record_id)

    def register_all_files(self, file_list):
        for fasta_file in file_list:
            for record in SeqIO.parse(fasta_file, "fasta"):
                self.register_record(record)

    def register_record(self, record):
        self.records.append(FancyRecord(record, self.tmpdir))

    def set_matrix(self, matrix):
        self.matrix_data = matrix
        for i in range(len(self.matrix_data)):
            record_i = self._get_record(self.matrix_data[i][0]["i"])
            for j in range(len(self.matrix_data[i])):
                record_j = self._get_record(self.matrix_data[i][j]["j"])
                self.matrix_data[i][j]["subplot"] = Subplot(
                    record_i, record_j, self.tmpdir, self.zoom
                )

        # More processing?
        logging.debug("\n" + pprint.pformat(matrix))

    def generate_matrix(self, mtype="complete"):
        matrix = []
        if mtype == "complete":
            for i in self.records:
                row = []
                for j in self.records:
                    row.append({"i": i.id, "j": j.id})
                matrix.append(row)
        elif mtype == "1vn":
            if len(self.records) < 2:
                raise RuntimeError("1vN not available for fewer than two sequences")
            else:
                row = []
                for j in self.records[1:]:
                    row.append({"i": self.records[0].id, "j": j.id})
                matrix.append(row)
        return matrix

    @classmethod
    def obtain_image_dimensions(cls, path):
        cmd = ["identify", path]
        output = subprocess.check_output(cmd, universal_newlines=True)
        size = output.split(" ")[3]
        (w, h) = size[0 : size.index("+")].split("x")
        return (int(w), int(h))

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
        cmd = ["convert", "-resize", scale, from_file, to_file]
        log.debug(" ".join(cmd))
        subprocess.check_call(cmd)

    def get_image_map(self):
        image_template = Template(
            '<area shape="rect" coords="${x1},${y1},${x2},${y2}" alt="${alt}" href="${href}" />'
        )
        imagemap = []

        j_widths = []
        i_height = []

        for j in range(len(self.matrix_data[0])):
            j_widths.append(self.matrix_data[0][j]["subplot"].get_thumb_dims()[0])

        for i in range(len(self.matrix_data)):
            i_height.append(self.matrix_data[i][0]["subplot"].get_thumb_dims()[1])

        log.debug(pprint.pformat(j_widths))
        log.debug(pprint.pformat(i_height))

        def cur_y(i_idx):
            return (
                MONTAGE_BORDER
                + sum(i_height[0:i_idx])
                + (2 * IMAGE_BORDER * (1 + i_idx))
            )

        def cur_x(j_idx):
            return (
                MONTAGE_BORDER
                + sum(j_widths[0:j_idx])
                + (2 * IMAGE_BORDER * (1 + j_idx))
            )

        for j in range(len(self.matrix_data[0])):
            for i in range(len(self.matrix_data)):
                current = self.matrix_data[i][j]["subplot"]
                # Move to final resting place
                new_image_location = current.move_to_dir(self.files_path)

                # Build imagemagp string
                imagemap.append(
                    image_template.substitute(
                        {
                            # Start at +image border so the border isn't included in
                            # start of box
                            "x1": cur_x(j),
                            "y1": cur_y(i),
                            "x2": cur_x(j) + j_widths[j],
                            "y2": cur_y(i) + i_height[i],
                            "alt": current.get_description(),
                            "href": new_image_location,
                        }
                    )
                )
        return "\n".join(imagemap)

    def _generate_montage(self):
        image_list = []
        for i in range(len(self.matrix_data)):
            for j in range(len(self.matrix_data[i])):
                subplot = self.matrix_data[i][j]["subplot"]
                image_list.append(subplot.thumb_path)

        # Montage step
        m0 = os.path.join(self.tmpdir, "m0.png")
        cmd = ["montage"] + image_list
        cmd += [
            "-tile",
            "%sx%s" % (len(self.matrix_data[0]), len(self.matrix_data)),
            "-geometry",
            "+0+0",
            "-border",
            str(IMAGE_BORDER),
            "-bordercolor",
            IMAGE_BORDER_COLOUR,
            "-font",
            TYPEFONT,
            m0,
        ]

        log.debug(" ".join(cmd))
        subprocess.check_call(cmd)

        # Add grey borders
        montage_path = os.path.join(self.tmpdir, "montage.png")
        cmd = [
            "convert",
            m0,
            "-bordercolor",
            IMAGE_BORDER_COLOUR,
            "-border",
            IMAGE_BORDER_COORD,
            "-bordercolor",
            MONTAGE_BORDER_COLOUR,
            "-border",
            MONTAGE_BORDER_COORD,
            montage_path,
        ]

        log.debug(" ".join(cmd))
        subprocess.check_call(cmd)
        os.unlink(m0)
        return montage_path

    def _annotate_montage(self, base_path):
        # Calculate some expected dimension
        cumulative_width = 0
        cumulative_height = 0
        for i in range(len(self.matrix_data)):
            for j in range(len(self.matrix_data[i])):
                subplot = self.matrix_data[i][j]["subplot"]

                if i == 0:
                    cumulative_width += subplot.get_thumb_dims()[0] + IMAGE_BORDER * 2

                if j == 0:
                    cumulative_height += subplot.get_thumb_dims()[1] + IMAGE_BORDER * 2

        convert_arguments_top = []
        convert_arguments_left = []
        left_offset = cumulative_width + MONTAGE_BORDER

        current_sum_width = MONTAGE_BORDER
        current_sum_height = MONTAGE_BORDER

        # Top side
        for j in range(len(self.matrix_data[0])):
            subplot = self.matrix_data[0][j]["subplot"]
            convert_arguments_top += [
                "-fill",
                LABEL_COLOUR,
                "-annotate",
                "+%s+40" % current_sum_width,
                subplot.j.header,
            ]
            current_sum_width += subplot.get_thumb_dims()[0] + (2 * IMAGE_BORDER)
            log.debug("CSW %s", current_sum_width)

        # Left side
        for i in range(len(self.matrix_data)):
            subplot = self.matrix_data[i][0]["subplot"]
            convert_arguments_left += [
                "-fill",
                LABEL_COLOUR,
                "-annotate",
                "+%s+%s" % (current_sum_height, left_offset),
                "\n" + subplot.i.header,
            ]
            current_sum_height += subplot.get_thumb_dims()[1] + (2 * IMAGE_BORDER)
            log.debug("CSH %s", current_sum_height)

        cmd = [
            "convert",
            base_path,
            "-rotate",
            "-90",
            "-pointsize",
            "20",
            "-font",
            TYPEFONT,
        ]
        cmd += convert_arguments_left
        cmd += ["-rotate", "90"]
        cmd += convert_arguments_top

        output_path = os.path.join(self.tmpdir, "large.png")
        cmd += [
            "-pointsize",
            "14",
            "-annotate",
            "+%s+%s" % (2, current_sum_height + 15),
            CREDITS,
            output_path,
        ]
        log.debug(" ".join(cmd))
        subprocess.check_output(cmd)
        return output_path

    def run(self):
        # We want the final thumbnail for the overall image to have a max width/height
        # of 700px for nice display
        total_seqlen_j = 0
        total_seqlen_i = 0
        for j in range(len(self.matrix_data[0])):
            total_seqlen_j += self.matrix_data[0][j]["subplot"].j.length

        for i in range(len(self.matrix_data)):
            total_seqlen_i += self.matrix_data[i][0]["subplot"].i.length
        total_seqlen = max(total_seqlen_i, total_seqlen_j)
        rescale = 200 * (700.0 / (float(total_seqlen) / self.zoom))
        rescale_p = "%s%%" % rescale

        # Generate gepard plots for each of the sub-images
        for i in range(len(self.matrix_data)):
            for j in range(len(self.matrix_data[i])):
                subplot = self.matrix_data[i][j]["subplot"]

                # generates _orig and _thumb versions
                subplot.run_gepard(self.matrix, self.window, global_rescale=rescale_p)

        base_montage = self._generate_montage()
        annotated_montage = self._annotate_montage(base_montage)

        final_montage_path = os.path.join(self.files_path, "large.png")
        shutil.move(annotated_montage, final_montage_path)
        return final_montage_path


def mist_wrapper(
    files,
    window=10,
    zoom=50,
    matrix="edna",
    plot_type="complete",
    files_path="mist_images",
):
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

    m = Misty(window=window, zoom=zoom, matrix=matrix, files_path=files_path)

    # There is only one special case, so we handle that separately. Every other
    # plot type wants ALL of the sequences available.
    if (
        plot_type == "2up"
        and len(files) != 2
        and matrix not in ("protidentity", "blosum62")
    ):
        idx = 0
        # Pull top two sequences.
        for fasta_file in files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                m.register_record(record)

                # Exit after we've seen two sequences
                idx += 1
                if idx == 2:
                    break
    else:
        m.register_all_files(files)

    # Generate the matrix appropriate to this plot type. There are two matrix
    # types: 1vN and complete. 1vN is just a line, complete is a complete
    # square.
    if plot_type == "complete":
        # ALL sequences are used.
        m.set_matrix(m.generate_matrix(mtype="complete"))
    elif plot_type == "2up":
        # Just two sequences.
        if len(files) == 2 and matrix in ("protidentity", "blosum62"):
            m.set_matrix(m.generate_matrix(mtype="complete"))
        else:
            # Co-op the 1vn to do a single plot, rather than a "complete" plot
            m.set_matrix(m.generate_matrix(mtype="1vn"))
    elif plot_type == "1vn":
        m.set_matrix(m.generate_matrix(mtype="1vn"))
    else:
        raise ValueError("Unknown plot type %s" % plot_type)

    m.run()
    # image_map will be returned from MIST
    image_map = m.get_image_map()

    return html_page % image_map


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MIST v3", epilog="")
    parser.add_argument(
        "files", nargs="+", type=argparse.FileType("r"), help="Fasta sequences"
    )

    parser.add_argument("--zoom", type=int, help="# bases / pixel", default=50)
    parser.add_argument("--window", type=int, help="Window size", default=10)
    parser.add_argument(
        "--matrix",
        type=str,
        choices=["ednaorig", "pam250", "edna", "protidentity", "blosum62"],
        help="# bases / pixel",
        default="edna",
    )

    parser.add_argument(
        "--plot_type",
        choices=["2up", "1vn", "complete"],
        help="Plot type",
        default="complete",
    )

    parser.add_argument(
        "--files_path", type=str, help="Image directory", default="mist_images"
    )

    args = parser.parse_args()
    print(mist_wrapper(**vars(args)))
