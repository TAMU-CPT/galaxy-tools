#!/usr/bin/env python
import argparse
from Bio import SeqIO
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name="charges")

HTML_HEADER = "<html><head><title>Charges Report</title></head><body>"
HTML_FOOTER = "</body></html>"

SVG_HEADER = '<svg width="%i" height="%i" xmlns="http://www.w3.org/2000/svg">\n'  # % (calcWidth, calcHeight)
SVG_FOOTER = "</svg>"

FULL_AA = [
    "H",
    "S",
    "Q",
    "T",
    "N",
    "C",
    "Y",
    "A",
    "V",
    "I",
    "L",
    "M",
    "P",
    "F",
    "W",
    "G",
    "E",
    "R",
    "D",
    "K",
]


def charges_html(svg, fasta, aa, fgColor, bgColor, width=120):
    colour_scheme = zip([x.upper() for x in aa], bgColor, fgColor)

    # CSS and header styling
    css = """<style type="text/css">
    .list li { list-style: none; margin:10px }
    .info { float:left; width:20px }
    pre { font-size:1.3em }
    """
    info = '<h1>Charges</h1><h3>Legend</h3><ul class="list">'
    for group in colour_scheme:
        css += ".%s{ background: %s; color: %s}\n" % group
        info += '<li><span class="%s" style="padding:5px">%s</span></li>\n' % (
            group[0],
            group[0],
        )
    css += "</style>"
    info += "</ul>"

    # Pre-calculate, so we can use for testing 'in'
    match_list = [group[0] for group in colour_scheme]

    page = ""
    # Parse sequences from fasta file
    for record in SeqIO.parse(fasta, "fasta"):
        page += "<pre><h3>&gt;%s %s</h3>\n" % (record.id, record.description)
        seq = list(str(record.seq).upper())

        idx = 0
        for i in range(0, len(seq), width):
            line_charges = []
            line_residues = seq[i : i + width]
            line_numbers = []

            for char in range(len(line_residues)):
                if line_residues[char] in "KRkr":
                    line_charges.append("+")
                elif line_residues[char] in "DEde":
                    line_charges.append("-")
                else:
                    line_charges.append(" ")

                # Could be swapped out for math with i+char...
                idx += 1
                if idx % 10 == 0:
                    line_numbers.append("%10s" % idx)

                # Replace with <span>
                for m in match_list:
                    if line_residues[char].upper() in m:
                        line_residues[char] = '<span class="%s">%s</span>' % (
                            m,
                            line_residues[char],
                        )

            page += "".join(line_charges) + "\n"
            page += "".join(line_residues) + "\n"
            page += "".join(line_numbers) + "\n"
            page += "\n"
        page += "</pre>"
    return HTML_HEADER + css + info + page + HTML_FOOTER


def charges_svg(svg, fasta, aa, fgColor, bgColor, width=120):
    colour_scheme = zip([x.upper() for x in aa], bgColor, fgColor)

    svgWidth = 1100

    # CSS and header styling
    classList = []
    classes = '<style type="text/css">\n<![CDATA[\n'

    defClass = ""
    for x in FULL_AA:
        addAA = True
        for y in aa:
            if x in y:
                addAA = False
        if addAA:
            defClass += x

    defBox = "#ffffff"
    defText = "#000000"

    for group in colour_scheme:
        classList.append(group[0])
        classes += "text.text_%s{fill: %s;}\n" % (group[0], group[2])
        classes += "rect.rect_%s{fill: %s; stroke: %s;}\n" % (
            group[0],
            group[1],
            group[1],
        )
        # info += '<li><span class="%s" style="padding:5px">%s</span></li>\n' % (group[0], group[0])
    if defClass != "":
        classes += "text.text_%s{fill: %s;}\n" % (defClass, defText)
        classes += "rect.rect_%s{fill: %s; stroke: %s;}\n" % (defClass, defBox, defBox)
        classList.append(defClass)
    classes += "text.info_text{white-space: pre;}\n"
    classes += "rect.rEven{fill: #fdfdfd; stroke: #fbfbfb;}\n"
    classes += "rect.rOdd{fill: #f2f2fc; stroke: #fbfbfb;}\n"
    classes += "]]></style>\n"
    body = ""
    groups = ""
    # Pre-calculate, so we can use for testing 'in'

    match_list = aa
    prevIndex = -1
    boxLen = 0
    page = ""
    title = ""

    yInd = 60
    yInc = 15
    seqIndent = 35
    idIndent = 20
    letterLen = 8.4375
    recNum = -1

    title += (
        '<text x="'
        + str(idIndent * 0.5)
        + '" y="'
        + str(yInd)
        + '" style="font-weight:bold; font-size:40px">Charges</text>\n'
    )
    yInd += 2 * yInc
    title += (
        '<text x="'
        + str(idIndent)
        + '" y="'
        + str(yInd)
        + '" style="font-size:18px">Legend:</text>\n'
    )
    yInd += 2 * yInc

    for i in range(len(classList)):
        title += (
            '<rect x="'
            + str(seqIndent)
            + '" y="'
            + str(yInd - yInc + 2)
            + '" width="'
            + str(len(classList[i]) * letterLen)
            + '" height="'
            + str(yInc)
            + '" class="rect_%s"/>\n' % classList[i]
        )
        title += (
            '<text x="'
            + str(seqIndent)
            + '" y="'
            + str(yInd)
            + '" class="text_%s" font-family="monospace" font-size="14">%s</text>\n'
            % (classList[i], classList[i])
        )
        yInd += yInc + 3
    yInd += yInc * 1.5

    # Parse sequences from fasta file
    for record in SeqIO.parse(fasta, "fasta"):

        recNum += 1
        seqHeader = (
            '<g><text x="'
            + str(idIndent)
            + '" y="'
            + str(yInd)
            + '" style="font-weight:bold">&gt;%s %s</text>\n'
            % (record.id, record.description)
        )
        body += seqHeader
        seq = list(str(record.seq).upper())
        yTop = yInd - yInc - 3
        yInd += yInc
        idx = 0
        for i in range(0, len(seq), width):
            line_charges = []
            line_residues = seq[i : i + width]
            line_numbers = []

            boxList = []
            groupList = []
            seqList = []
            prevIndex = -1
            boxLen = 0
            for char in range(len(line_residues)):

                thisInd = 0
                for x in match_list:

                    if line_residues[char] in x:
                        break
                    thisInd += 1

                if thisInd == len(match_list):
                    thisInd = -1

                if char != 0 and thisInd != prevIndex:
                    boxList.append(boxLen)
                    seqList.append((line_residues[char - boxLen : char]))
                    groupList.append(prevIndex)
                    boxLen = 0
                prevIndex = thisInd
                boxLen += 1

                if line_residues[char] in "KRkr":
                    line_charges.append("+")
                elif line_residues[char] in "DEde":
                    line_charges.append("-")
                else:
                    line_charges.append(" ")

                # Could be swapped out for math with i+char...
                idx += 1
                if idx % 10 == 0:
                    line_numbers.append("%10s" % idx)

                # Replace with <span>
                # for m in match_list:
                #    if line_residues[char].upper() in m:
                #        line_residues[char] = '<span class="%s">%s</span>' % (m, line_residues[char])

            seqList.append((line_residues[-boxLen:]))
            boxList.append(boxLen)
            groupList.append(prevIndex)
            # Write line charges
            line = "".join(line_charges)
            body += (
                '<text x="'
                + str(seqIndent)
                + '" y="'
                + str(yInd)
                + '" class="info_text" font-family="monospace" font-size="14">%s</text>\n'
                % line
            )
            yInd += yInc
            # Write sequence
            sumSeq = 0
            for i in range(len(seqList)):
                res = ""
                for sub in seqList[i]:
                    res += sub
                body += (
                    '<rect x="'
                    + str(.5 + seqIndent + (letterLen * sumSeq))
                    + '" y="'
                    + str(yInd - yInc + 2)
                    + '" width="'
                    + str(boxList[i] * letterLen)
                    + '" height="'
                    + str(yInc)
                    + '" class="rect_%s"/>\n' % classList[groupList[i]]
                )
                body += (
                    '<text x="'
                    + str(seqIndent + (letterLen * sumSeq))
                    + '" y="'
                    + str(yInd)
                    + '" class="text_%s" font-family="monospace" font-size="14">%s</text>\n'
                    % (classList[groupList[i]], res)
                )
                sumSeq += len(seqList[i])
            yInd += yInc
            # Write numbers
            line = "".join(line_numbers) + "\n"
            body += (
                '<text x="'
                + str(seqIndent)
                + '" y="'
                + str(yInd)
                + '" class="info_text" font-size="14" font-family="monospace">%s</text>\n'
                % line
            )
            yInd += yInc

        body += "</g>\n"
        yInd += yInc
        if recNum % 2 == 0:
            groups += (
                '<rect x="0" y="'
                + str(yTop)
                + '" width="'
                + str(svgWidth + 1)
                + '" height="'
                + str(yInd - yTop)
                + '" class="rEven"/>\n'
            )
        else:
            groups += (
                '<rect x="0" y="'
                + str(yTop)
                + '" width="'
                + str(svgWidth + 1)
                + '" height="'
                + str(yInd - yTop)
                + '" class="rOdd"/>\n'
            )
    svgHeight = yInd

    return (
        (SVG_HEADER % (svgWidth, svgHeight))
        + title
        + classes
        + groups
        + body
        + SVG_FOOTER
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Top related genomes")
    parser.add_argument("--svg", action="store_true")
    parser.add_argument("fasta", type=argparse.FileType("r"), help="Fasta protein file")
    parser.add_argument("--width", type=int, help="Plot width", default=120)
    parser.add_argument("--aa", nargs="+")
    parser.add_argument("--fgColor", nargs="+")
    parser.add_argument("--bgColor", nargs="+")

    args = parser.parse_args()
    if args.svg:
        print(charges_svg(**vars(args)))
    else:
        print(charges_html(**vars(args)))
