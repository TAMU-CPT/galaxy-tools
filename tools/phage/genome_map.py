#!/usr/bin/env python
# vim: set fileencoding=utf-8
import base64
import argparse
import svgwrite
import logging
from gff3 import feature_lambda, feature_test_type, get_gff3_id, wa_unified_product_name
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, ExactPosition
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

FONT_STYLE = "fill: #000033; stroke:none;font-family: 'monospace';font-size:9px;"
DEFAULT_COLOR_SCHEME = {
    "tRNA": {
        "color": "#ee0000",
        "border": 0,
        "plot": 1
    },
    "gene": {
        "color": "#0086ee",
        "border": 1,
        "plot": 1
    },
    "CDS": {
        "color": "#1c86ee",
        "border": 1,
        "plot": 1
    },
    "mRNA": {
        "color": "#ff0000",
        "border": 1,
        "plot": 1
    },
    "regulatory": {
        "color": "#000000",
        "border": 0,
        "plot": 1
    },
    "repeat_region": {
        "color": "#b3ee3a",
        "border": 1,
        "plot": 1
    },
    "mat_peptide": {
        "color": "#b23aee",
        "border": 1,
        "plot": 1
    }
}


class PlottedFeature(object):

    def __init__(self, feature):
        self.feature = feature
        self.location = feature.location
        self.tag = feature.type

        def featureColor(feature, allow_default=True):
            if 'color' in feature.qualifiers:
                color = feature.qualifiers['color'][0]
            elif 'colour' in feature.qualifiers:
                color = feature.qualifiers['colour'][0]
            elif allow_default and feature.type in DEFAULT_COLOR_SCHEME:
                color = DEFAULT_COLOR_SCHEME[feature.type]['color']
            else:
                color = None
            return color

        non_default_cds_colors = [featureColor(x) for x in self.get_cdss() if featureColor(x, allow_default=False)]
        if len(non_default_cds_colors) > 0:
            self.color = non_default_cds_colors[0]
        else:
            self.color = featureColor(feature, allow_default=True) or '#000000'

    def get_cdss(self):
        return list(feature_lambda(self.feature.sub_features, feature_test_type, {'type': 'CDS'}, subfeatures=False))

    def get_label(self):
        if hasattr(self, 'label'):
            return self.label

        label = wa_unified_product_name(self.feature)
        # Best option, the gene has a product
        if 'product' in self.feature.qualifiers or 'Product' in self.feature.qualifiers:
            self.label = label
            return self.label

        # Otherwise, fail over the the CDSs
        cdss = self.get_cdss()
        labels = [wa_unified_product_name(cds) for cds in cdss]

        # If the CDSs have a label, use that.
        if len(labels) > 0:
            self.label = labels[0]
        else:
            # Or fail over to empty
            self.label = None

        return self.label


class GeneClass(object):

    def __init__(self, objects=None, key='gene', color='blue', border=1, plot=True, included=True):
        self.objects = objects
        if objects is None:
            self.objects = []
        self.included = included
        self.key = key
        self.color = color
        self.border = border
        self.plot = plot

    def addObject(self, obj):
        self.objects += [obj]

    def getMemberCount(self):
        return len(self.objects)


class Plotter(object):

    def __init__(self, rows=2, hypo=False):
        self.line_Count = 1
        self._ft_count = 0

        self.justified = True
        self.separate_strands = True
        self.double_line_for_overlap = True
        self.opacity = 1.0

        self.label_hypo = hypo

        self.genome_length = 0

        self.x_offset = 30
        self.y_offset = 90
        self.ils = 200
        self.zoom = 80
        self.rows = rows

        self.rowdata = []
        self._internal_maxrowlength = 0

        self.cs = DEFAULT_COLOR_SCHEME
        self.classes = {}
        for key in self.cs:
            self.classes[key] = GeneClass(
                key=key,
                color=DEFAULT_COLOR_SCHEME[key]['color'],
                border=DEFAULT_COLOR_SCHEME[key]['border'],
                plot=DEFAULT_COLOR_SCHEME[key]['plot'],
                included=True,
            )

    def processSequence(self, seq):
        self.genome_length = len(seq)

        if self.genome_length < 50000:
            log.warn("Sequence shorter than 10kb, capped at 1 row")
            self.rows = 1

        for feature in seq.features:
            self.processFeature(feature)

        self.ruler_offset = seq.subset

    def processFeature(self, feature):
        gene = PlottedFeature(feature)

        if feature.type in self.classes:
            self.classes[feature.type].addObject(gene)

    def partitionLines(self, split_factor=1.05):
        avgRowLength = int(float(self.genome_length) / float(self.rows * split_factor))

        fake_count = 100

        items = []
        for i in range(fake_count):
            key = int(float(self.genome_length * i) / fake_count)
            items.append(FeatureLocation(key, key, strand=1))

        for x in self.classes:
            if self.classes[x].included:
                items += [y.location for y in self.classes[x].objects]

        longest_last_object = 1
        thisRowEnd = 1 + avgRowLength
        currentRow = 1
        _internal_maxrowlength = 0

        rowData = {
            1: {
                'start': ExactPosition(1)
            }
        }

        for item in sorted(items, key=lambda x: x.start):
            if item.start >= thisRowEnd or item.end > thisRowEnd:

                if self.justified or item.start >= rowData[currentRow]['end']:
                    rowData[currentRow]['end'] = thisRowEnd
                else:
                    rowData[currentRow]['end'] = max(longest_last_object, item.start)

                _internal_maxrowlength = max(
                    _internal_maxrowlength,
                    rowData[currentRow]['end'] - rowData[currentRow]['start']
                )

                currentRow += 1
                rowData[currentRow] = {}

                if item.start <= rowData[currentRow - 1]['end']:
                    rowData[currentRow]['start'] = item.start
                else:
                    rowData[currentRow]['start'] = rowData[currentRow - 1]['end'] + 1

                thisRowEnd = avgRowLength + rowData[currentRow]['start']

        thisRowEnd = rowData[currentRow]['end'] = ExactPosition(self.genome_length + 1)

        _internal_maxrowlength = max(
            _internal_maxrowlength,
            rowData[currentRow]['end'] - rowData[currentRow]['start']
        )

        return rowData, avgRowLength, _internal_maxrowlength

    def createSvg(self, rowData, widthOverride=0):
        height = int((len(rowData.keys())) * self.ils)
        width = int(float(self.avgRowLength) / self.zoom)

        if widthOverride != 0:
            width = widthOverride

        self.calc_width = width

        self.svg = svgwrite.Drawing(
            size=(
                width + 2 * (self.x_offset),
                height + 2 * (self.y_offset)
            )
        )

        ui_group = self.svg.g(
            id="group_ui",
            style="stroke: #000000; fill: #000000; fill-opacity: 1;"
        )

        for i in range(1, max(rowData.keys()) + 1):
            self.addRuler(i, ui_group, rowData)

        for key in self.classes:
            c = self.classes[key]

            if not c.plot:
                continue

            class_group = self.svg.g(
                id="group_%s" % key,
                style="stroke: %s; fill: %s; fill-opacity: %s" % (
                    'black' if c.border else 'none', c.color, self.opacity)
            )

            for gene in c.objects:
                svgFeature, x, y, w, h = self.featureBox(
                    gene,
                    rowData,
                    class_group,
                    self.calculateRow(gene, rowData),
                )
                class_group.add(svgFeature)

                if gene.get_label():
                    if self.label_hypo or (not self.label_hypo and 'ypothetical' not in gene.get_label()):
                        # print self.label_hypo, (self.label_hypo and 'ypothetical' not in gene.get_label()), 'ypothetical' not in gene.get_label(), gene.get_label()
                        # if not self.label_hypo or (self.label_hypo and 'ypothetical' not in gene.get_label()):
                        svgFeatureLabel = self.featureLabel(
                            gene, rowData, class_group,
                            self.calculateRow(gene, rowData),
                            gene.get_label(),
                            x, y, w, h
                        )
                        class_group.add(svgFeatureLabel)

            self.svg.add(class_group)

        self.svg.add(ui_group)
        return self.svg

    def fitness(self, rowData, maxRowLength):
        score = 2
        score -= abs(len(rowData.keys()) - self.rows)

        lastRow = rowData[max(rowData.keys())]
        lastRow = float(lastRow['end'] - lastRow['start'])
        lastRow /= maxRowLength

        return score + lastRow

    def optimizedPartition(self):
        bestRowData = None
        bestFitness = 0
        for i in range(70, 200):
            s = float(i) / 100
            results = self.partitionLines(split_factor=s)
            fitness = self.fitness(results[0], results[2])
            # print s, fitness
            if fitness > bestFitness:
                bestRowData = results
                bestFitness = fitness

        self.avgRowLength = bestRowData[1]
        self._internal_maxrowlength = bestRowData[2]
        return bestRowData[0]

    def offsetPoint(self, x, y):
        return (x + self.x_offset, y + self.y_offset)

    def addRuler(self, row, ui_group, rowData):
        y_fix = self.ils * (row - 1)

        line_width = self.calc_width * (rowData[row]['end'] - rowData[row]['start']) / self._internal_maxrowlength

        ui_group.add(self.svg.line(
            start=self.offsetPoint(0, y_fix),
            end=self.offsetPoint(line_width, y_fix),
            id='ruler_%s' % row,
        ))

        if self.double_line_for_overlap and row > 1:
            if rowData[row - 1]['end'] - rowData[row]['start'] >= 0:
                line_length = self.calc_width * (rowData[row - 1]['end'] - rowData[row]['start']) / self._internal_maxrowlength

                if line_length > 0:

                    ui_group.add(self.svg.line(
                        id='ruler_%s_overlap' % row,
                        start=self.offsetPoint(
                            0,
                            y_fix - 5
                        ),
                        end=self.offsetPoint(
                            line_length,
                            y_fix - 5
                        )
                    ))

        # ui_group.add(svg.line(
            # id='ruler_%s_left_border' % row,
            # start=self.offsetPoint(0, y_fix),
            # end = self.offsetPoint(line_width, y_fix)
        # ))

        for idx in range(rowData[row]['start'] - 1, rowData[row]['end']):
            if idx % 1000 == 0:
                current_location = self.calc_width * (idx - rowData[row]['start']) / self._internal_maxrowlength

                line_height = 5
                if idx % 10000 == 0:
                    line_height = 10

                ui_group.add(self.svg.line(
                    id='ruler_vert_%s' % idx,
                    start=self.offsetPoint(current_location, y_fix),
                    end=self.offsetPoint(current_location, y_fix + line_height),
                ))

                if idx % 10000 == 0:
                    ui_group.add(self.svg.text(
                        id='ruler_text_%s' % idx,
                        text='%s kb' % int((idx + self.ruler_offset) / 1000),
                        x=[current_location + 10 + self.x_offset],
                        y=[y_fix + 20 + self.y_offset],
                        style=FONT_STYLE,
                    ))

    def calculateRow(self, obj, rowData):
        for i in rowData.keys():
            if rowData[i]['start'] - 1 <= obj.location.start < rowData[i]['end'] and \
                    rowData[i]['start'] - 1 <= obj.location.end < rowData[i]['end']:
                return i
        raise Exception("Cannot place feature")

    def featureBox(self, feature, rowData, class_group, row):
        x = self.calc_width * (
            feature.location.start - rowData[row]['start']
        ) / self._internal_maxrowlength + self.x_offset
        h = 15
        y = (row - 1) * self.ils + self.y_offset - h / 2

        w = self.calc_width * len(feature.location) / self._internal_maxrowlength

        if self.separate_strands:
            y += -30 * feature.location.strand

        # Alternate acording to frame
        y += \
            10 * ((feature.location.start - 2 * feature.location.strand + 1) % 3) - \
            10 * feature.location.strand

        return self.svg.rect(
            id=base64.b32encode(get_gff3_id(feature.feature)),
            insert=(x, y),
            size=(w, h),
            style="fill:%s;stroke-width:0.5;" % feature.color,
        ), x, y, w, h

    def featureLabel(self, feature, rowData, class_group, row, label, x, y, w, h):
        lx = x + w / 2 - len(label) * 4
        lxm = x + (float(w) / 2)
        ly = y + h / 2 + 40

        if feature.location.strand > 0:
            ly -= 70

        g = self.svg.g(
            id="%s_g" % base64.b32encode(get_gff3_id(feature.feature)),
        )
        g.add(self.svg.text(
            id='label_text_%s' % base64.b32encode(get_gff3_id(feature.feature)),
            text=label,
            x=[lx],
            y=[ly],
            style=FONT_STYLE,
        ))

        if feature.location.strand > 0:
            callout_start = (lxm, y)
            callout_end = (lxm, ly)
        else:
            callout_start = (lxm, y + h)
            callout_end = (lxm, ly - h)

        g.add(self.svg.line(
            id='label_callout_%s' % base64.b32encode(get_gff3_id(feature.feature)),
            start=callout_start,
            end=callout_end,
        ))

        return g


def parseFile(annotations, genome, subset=None, rows=2, width=0, hypo=False):
    plotter = Plotter(rows=rows, hypo=hypo)

    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    for record in GFF.parse(annotations, base_dict=seq_dict):
        if subset is not None:
            (a, b) = map(int, subset.split(','))
            record = record[a:b]
            record.subset = a
        else:
            record.subset = 0

        plotter.processSequence(record)
        rowData = plotter.optimizedPartition()
        svg = plotter.createSvg(rowData, widthOverride=width)
        print svg.tostring()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations')
    parser.add_argument('annotations', type=file, help='Parent GFF3 annotations')
    parser.add_argument('genome', type=file, help='Genome Sequence')
    parser.add_argument('--subset', help="Subset location (E.g. --subset '100,400')")
    parser.add_argument('--rows', default=2, type=int, help="Number of rows")
    parser.add_argument('--width', default=0, type=int, help='Width of plot')
    parser.add_argument('--hypo', action='store_true', help='Label hypotheticals')
    args = parser.parse_args()

    parseFile(**vars(args))
