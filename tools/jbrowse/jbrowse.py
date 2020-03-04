#!/usr/bin/env python
import os
import copy
import argparse
import subprocess
import hashlib
import binascii
import struct
import datetime
import tempfile
import shutil
import json
from Bio.Data import CodonTable
import xml.etree.ElementTree as ET
import logging
from collections import defaultdict

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("jbrowse")
TODAY = datetime.datetime.now().strftime("%Y-%m-%d")
GALAXY_INFRASTRUCTURE_URL = None


class ColorScaling(object):

    COLOR_FUNCTION_TEMPLATE = """
    function(feature, variableName, glyphObject, track) {{
        var score = {score};
        {opacity}
        return 'rgba({red}, {green}, {blue}, ' + opacity + ')';
    }}
    """

    COLOR_FUNCTION_TEMPLATE_QUAL = """
    function(feature, variableName, glyphObject, track) {{
        var search_up = function self(sf, attr){{
            if(sf.get(attr) !== undefined){{
                return sf.get(attr);
            }}
            if(sf.parent() === undefined) {{
                return;
            }}else{{
                return self(sf.parent(), attr);
            }}
        }};

        var search_down = function self(sf, attr){{
            if(sf.get(attr) !== undefined){{
                return sf.get(attr);
            }}
            if(sf.children() === undefined) {{
                return;
            }}else{{
                var kids = sf.children();
                for(var child_idx in kids){{
                    var x = self(kids[child_idx], attr);
                    if(x !== undefined){{
                        return x;
                    }}
                }}
                return;
            }}
        }};

        var color = ({user_spec_color} || search_up(feature, 'color') || search_down(feature, 'color') || {auto_gen_color});
        var score = (search_up(feature, 'score') || search_down(feature, 'score'));
        {opacity}
        if(score === undefined){{ opacity = 1; }}
        var result = /^#?([a-f\d]{{2}})([a-f\d]{{2}})([a-f\d]{{2}})$/i.exec(color);
        var red = parseInt(result[1], 16);
        var green = parseInt(result[2], 16);
        var blue = parseInt(result[3], 16);
        if(isNaN(opacity) || opacity < 0){{ opacity = 0; }}
        return 'rgba(' + red + ',' + green + ',' + blue + ',' + opacity + ')';
    }}
    """

    OPACITY_MATH = {
        "linear": """
            var opacity = (score - ({min})) / (({max}) - ({min}));
        """,
        "logarithmic": """
            var opacity = (score - ({min})) / (({max}) - ({min}));
            opacity = Math.log10(opacity) + Math.log10({max});
        """,
        "blast": """
            var opacity = 0;
            if(score == 0.0) {{
                opacity = 1;
            }} else {{
                opacity = (20 - Math.log10(score)) / 180;
            }}
        """,
    }

    BREWER_COLOUR_IDX = 0
    BREWER_COLOUR_SCHEMES = [
        (166, 206, 227),
        (31, 120, 180),
        (178, 223, 138),
        (51, 160, 44),
        (251, 154, 153),
        (227, 26, 28),
        (253, 191, 111),
        (255, 127, 0),
        (202, 178, 214),
        (106, 61, 154),
        (255, 255, 153),
        (177, 89, 40),
        (228, 26, 28),
        (55, 126, 184),
        (77, 175, 74),
        (152, 78, 163),
        (255, 127, 0),
    ]

    BREWER_DIVERGING_PALLETES = {
        "BrBg": ("#543005", "#003c30"),
        "PiYg": ("#8e0152", "#276419"),
        "PRGn": ("#40004b", "#00441b"),
        "PuOr": ("#7f3b08", "#2d004b"),
        "RdBu": ("#67001f", "#053061"),
        "RdGy": ("#67001f", "#1a1a1a"),
        "RdYlBu": ("#a50026", "#313695"),
        "RdYlGn": ("#a50026", "#006837"),
        "Spectral": ("#9e0142", "#5e4fa2"),
    }

    def __init__(self):
        self.brewer_colour_idx = 0

    def rgb_from_hex(self, hexstr):
        # http://stackoverflow.com/questions/4296249/how-do-i-convert-a-hex-triplet-to-an-rgb-tuple-and-back
        return struct.unpack("BBB", binascii.unhexlify(hexstr))

    def min_max_gff(self, gff_file):
        min_val = None
        max_val = None
        with open(gff_file, "r") as handle:
            for line in handle:
                try:
                    value = float(line.split("\t")[5])
                    min_val = min(value, (min_val or value))
                    max_val = max(value, (max_val or value))

                    if value < min_val:
                        min_val = value

                    if value > max_val:
                        max_val = value
                except Exception:
                    pass
        return min_val, max_val

    def hex_from_rgb(self, r, g, b):
        return "#%02x%02x%02x" % (r, g, b)

    def _get_colours(self):
        r, g, b = self.BREWER_COLOUR_SCHEMES[
            self.brewer_colour_idx % len(self.BREWER_COLOUR_SCHEMES)
        ]
        self.brewer_colour_idx += 1
        return r, g, b

    def parse_colours(self, track, trackFormat, gff3=None):
        # Wiggle tracks have a bicolor pallete
        trackConfig = {"style": {}}
        if trackFormat == "wiggle":

            trackConfig["style"]["pos_color"] = track["wiggle"]["color_pos"]
            trackConfig["style"]["neg_color"] = track["wiggle"]["color_neg"]

            if trackConfig["style"]["pos_color"] == "__auto__":
                trackConfig["style"]["neg_color"] = self.hex_from_rgb(
                    *self._get_colours()
                )
                trackConfig["style"]["pos_color"] = self.hex_from_rgb(
                    *self._get_colours()
                )

            # Wiggle tracks can change colour at a specified place
            bc_pivot = track["wiggle"]["bicolor_pivot"]
            if bc_pivot not in ("mean", "zero"):
                # The values are either one of those two strings
                # or a number
                bc_pivot = float(bc_pivot)
            trackConfig["bicolor_pivot"] = bc_pivot
        elif "scaling" in track:
            if track["scaling"]["method"] == "ignore":
                if track["scaling"]["scheme"]["color"] != "__auto__":
                    trackConfig["style"]["color"] = track["scaling"]["scheme"]["color"]
                else:
                    trackConfig["style"]["color"] = self.hex_from_rgb(
                        *self._get_colours()
                    )
            else:
                # Scored method
                algo = track["scaling"]["algo"]
                # linear, logarithmic, blast
                scales = track["scaling"]["scales"]
                # type __auto__, manual (min, max)
                scheme = track["scaling"]["scheme"]
                # scheme -> (type (opacity), color)
                # ==================================
                # GENE CALLS OR BLAST
                # ==================================
                if trackFormat == "blast":
                    red, green, blue = self._get_colours()
                    color_function = self.COLOR_FUNCTION_TEMPLATE.format(
                        **{
                            "score": "feature._parent.get('score')",
                            "opacity": self.OPACITY_MATH["blast"],
                            "red": red,
                            "green": green,
                            "blue": blue,
                        }
                    )
                    trackConfig["style"]["color"] = color_function.replace("\n", "")
                elif trackFormat == "gene_calls":
                    # Default values, based on GFF3 spec
                    min_val = 0
                    max_val = 1000
                    # Get min/max and build a scoring function since JBrowse doesn't
                    if scales["type"] == "automatic" or scales["type"] == "__auto__":
                        min_val, max_val = self.min_max_gff(gff3)
                    else:
                        min_val = scales.get("min", 0)
                        max_val = scales.get("max", 1000)

                    if scheme["color"] == "__auto__":
                        user_color = "undefined"
                        auto_color = "'%s'" % self.hex_from_rgb(*self._get_colours())
                    elif scheme["color"].startswith("#"):
                        user_color = "'%s'" % self.hex_from_rgb(
                            *self.rgb_from_hex(scheme["color"][1:])
                        )
                        auto_color = "undefined"
                    else:
                        user_color = "undefined"
                        auto_color = "'%s'" % self.hex_from_rgb(*self._get_colours())

                    color_function = self.COLOR_FUNCTION_TEMPLATE_QUAL.format(
                        **{
                            "opacity": self.OPACITY_MATH[algo].format(
                                **{"max": max_val, "min": min_val}
                            ),
                            "user_spec_color": user_color,
                            "auto_gen_color": auto_color,
                        }
                    )

                    trackConfig["style"]["color"] = color_function.replace("\n", "")
        return trackConfig


def etree_to_dict(t):
    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(("@" + k, v) for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]["#text"] = text
        else:
            d[t.tag] = text
    return d


# score comes from feature._parent.get('score') or feature.get('score')

INSTALLED_TO = os.path.dirname(os.path.realpath(__file__))


def metadata_from_node(node):
    metadata = {}
    try:
        if len(node.findall("dataset")) != 1:
            # exit early
            return metadata
    except:
        return {}

    for (key, value) in node.findall("dataset")[0].attrib.items():
        metadata["dataset_%s" % key] = value

    for (key, value) in node.findall("history")[0].attrib.items():
        metadata["history_%s" % key] = value

    for (key, value) in node.findall("metadata")[0].attrib.items():
        metadata["metadata_%s" % key] = value

    for (key, value) in node.findall("tool")[0].attrib.items():
        metadata["tool_%s" % key] = value

    # Additional Mappings applied:
    metadata[
        "dataset_edam_format"
    ] = '<a target="_blank" href="http://edamontology.org/{0}">{1}</a>'.format(
        metadata["dataset_edam_format"], metadata["dataset_file_ext"]
    )
    metadata["history_user_email"] = '<a href="mailto:{0}">{0}</a>'.format(
        metadata["history_user_email"]
    )
    metadata[
        "history_display_name"
    ] = '<a target="_blank" href="{galaxy}/history/view/{encoded_hist_id}">{hist_name}</a>'.format(
        galaxy=GALAXY_INFRASTRUCTURE_URL,
        encoded_hist_id=metadata["history_id"],
        hist_name=metadata["history_display_name"],
    )
    metadata[
        "tool_tool"
    ] = '<a target="_blank" href="{galaxy}/datasets/{encoded_id}/show_params">{tool_id}</a>'.format(
        galaxy=GALAXY_INFRASTRUCTURE_URL,
        encoded_id=metadata["dataset_id"],
        tool_id=metadata["tool_tool_id"],
        tool_version=metadata["tool_tool_version"],
    )
    return metadata


class JbrowseConnector(object):
    def __init__(self, jbrowse, outdir, genomes, standalone=False, gencode=1):
        self.TN_TABLE = {
            "gff3": "--gff",
            "gff": "--gff",
            "bed": "--bed",
            "genbank": "--gbk",
        }

        self.cs = ColorScaling()
        self.jbrowse = jbrowse
        self.outdir = outdir
        self.genome_paths = genomes
        self.standalone = standalone
        self.gencode = gencode

        if standalone:
            self.clone_jbrowse(self.jbrowse, self.outdir)
        else:
            try:
                os.makedirs(self.outdir)
            except OSError:
                # Ignore if the folder exists
                pass

            try:
                os.makedirs(os.path.join(self.outdir, "data", "raw"))
            except OSError:
                # Ignore if the folder exists
                pass

        self.process_genomes()
        self.update_gencode()

    def update_gencode(self):
        table = CodonTable.unambiguous_dna_by_id[int(self.gencode)]
        trackList = os.path.join(self.outdir, "data", "trackList.json")
        with open(trackList, "r") as handle:
            trackListData = json.load(handle)

        trackListData["tracks"][0].update(
            {
                "codonStarts": table.start_codons,
                "codonStops": table.stop_codons,
                "codonTable": table.forward_table,
            }
        )

        with open(trackList, "w") as handle:
            json.dump(trackListData, handle, indent=2)

    def subprocess_check_call(self, command):
        log.debug("cd %s && %s", self.outdir, " ".join(command))
        subprocess.check_call(command, cwd=self.outdir)

    def _jbrowse_bin(self, command):
        return os.path.realpath(os.path.join(self.jbrowse, "bin", command))

    def process_genomes(self):
        for genome_node in self.genome_paths:
            # TODO: Waiting on https://github.com/GMOD/jbrowse/pull/884
            self.subprocess_check_call(
                [
                    "perl",
                    self._jbrowse_bin("prepare-refseqs.pl"),
                    "--fasta",
                    genome_node["path"],
                ]
            )

        # Generate name
        # self.subprocess_check_call([
        # 'perl', self._jbrowse_bin('generate-names.pl'),
        # '--hashBits', '16'
        # ])

    def _add_json(self, json_data):

        cmd = [
            "perl",
            self._jbrowse_bin("add-json.pl"),
            json.dumps(json_data),
            os.path.join("data", "trackList.json"),
        ]
        self.subprocess_check_call(cmd)

    def _add_track_json(self, json_data):
        if len(json_data.keys()) == 0:
            return

        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.write(json.dumps(json_data))
        tmp.close()
        cmd = [
            "perl",
            self._jbrowse_bin("add-track-json.pl"),
            tmp.name,
            os.path.join("data", "trackList.json"),
        ]
        self.subprocess_check_call(cmd)
        os.unlink(tmp.name)

    def _blastxml_to_gff3(self, xml, min_gap=10):
        gff3_unrebased = tempfile.NamedTemporaryFile(delete=False)
        cmd = [
            "python",
            os.path.join(INSTALLED_TO, "blastxml_to_gapped_gff3.py"),
            "--trim",
            "--trim_end",
            "--min_gap",
            str(min_gap),
            xml,
        ]
        log.debug("cd %s && %s > %s", self.outdir, " ".join(cmd), gff3_unrebased.name)
        subprocess.check_call(cmd, cwd=self.outdir, stdout=gff3_unrebased)
        gff3_unrebased.close()
        return gff3_unrebased.name

    def add_blastxml(self, data, trackData, blastOpts, **kwargs):
        gff3 = self._blastxml_to_gff3(data, min_gap=blastOpts["min_gap"])

        if "parent" in blastOpts and blastOpts["parent"] != "None":
            gff3_rebased = tempfile.NamedTemporaryFile(delete=False)
            cmd = ["python", os.path.join(INSTALLED_TO, "gff3_rebase.py")]
            if blastOpts.get("protein", "false") == "true":
                cmd.append("--protein2dna")
            cmd.extend([os.path.realpath(blastOpts["parent"]), gff3])
            log.debug("cd %s && %s > %s", self.outdir, " ".join(cmd), gff3_rebased.name)
            subprocess.check_call(cmd, cwd=self.outdir, stdout=gff3_rebased)
            gff3_rebased.close()

            # Replace original gff3 file
            shutil.copy(gff3_rebased.name, gff3)
            os.unlink(gff3_rebased.name)

        config = {
            "glyph": "JBrowse/View/FeatureGlyph/Segments",
            "category": trackData["category"],
        }

        clientConfig = trackData["style"]

        cmd = [
            "perl",
            self._jbrowse_bin("flatfile-to-json.pl"),
            "--gff",
            gff3,
            "--trackLabel",
            trackData["label"],
            "--key",
            trackData["key"],
            "--clientConfig",
            json.dumps(clientConfig),
            "--config",
            json.dumps(config),
            "--trackType",
            "BlastView/View/Track/CanvasFeatures",
        ]

        self.subprocess_check_call(cmd)
        os.unlink(gff3)

    def add_bigwig(self, data, trackData, wiggleOpts, **kwargs):
        dest = os.path.join("data", "raw", trackData["label"] + ".bw")
        cmd = ["cp", data, dest]
        self.subprocess_check_call(cmd)

        trackData.update(
            {
                "urlTemplate": os.path.join("..", dest),
                "storeClass": "JBrowse/Store/SeqFeature/BigWig",
                "type": "JBrowse/View/Track/Wiggle/Density",
            }
        )

        trackData["type"] = wiggleOpts["type"]
        trackData["variance_band"] = (
            True if wiggleOpts["variance_band"] == "true" else False
        )

        if "min" in wiggleOpts and "max" in wiggleOpts:
            trackData["min_score"] = wiggleOpts["min"]
            trackData["max_score"] = wiggleOpts["max"]
        else:
            trackData["autoscale"] = wiggleOpts.get("autoscale", "local")

        trackData["scale"] = wiggleOpts["scale"]

        self._add_track_json(trackData)

    def add_bam(self, data, trackData, bamOpts, bam_index=None, **kwargs):
        dest = os.path.join("data", "raw", trackData["label"] + ".bam")
        cmd = ["cp", "-s", os.path.realpath(data), dest]
        self.subprocess_check_call(cmd)

        cmd = ["cp", "-s", os.path.realpath(bam_index), dest + ".bai"]
        self.subprocess_check_call(cmd)

        trackData.update(
            {
                "urlTemplate": os.path.join("..", dest),
                "type": "JBrowse/View/Track/Alignments2",
                "storeClass": "JBrowse/Store/SeqFeature/BAM",
            }
        )

        self._add_track_json(trackData)

        if bamOpts.get("auto_snp", "false") == "true":
            trackData2 = copy.copy(trackData)
            trackData2.update(
                {
                    "type": "JBrowse/View/Track/SNPCoverage",
                    "key": trackData["key"] + " - SNPs/Coverage",
                    "label": trackData["label"] + "_autosnp",
                }
            )
            self._add_track_json(trackData2)

    def add_vcf(self, data, trackData, vcfOpts={}, **kwargs):
        dest = os.path.join("data", "raw", trackData["label"] + ".vcf")
        # ln?
        cmd = ["cp", "-s", data, dest]
        self.subprocess_check_call(cmd)
        cmd = ["bgzip", dest]
        self.subprocess_check_call(cmd)
        cmd = ["tabix", "-p", "vcf", dest + ".gz"]
        self.subprocess_check_call(cmd)

        trackData.update(
            {
                "urlTemplate": os.path.join("..", dest + ".gz"),
                "type": "JBrowse/View/Track/HTMLVariants",
                "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
            }
        )
        self._add_track_json(trackData)

    def add_features(self, data, format, trackData, gffOpts, metadata=None, **kwargs):
        cmd = [
            "perl",
            self._jbrowse_bin("flatfile-to-json.pl"),
            self.TN_TABLE.get(format, "gff"),
            data,
            "--trackLabel",
            trackData["label"],
            # '--trackType', 'JBrowse/View/Track/CanvasFeatures',
            "--key",
            trackData["key"],
        ]

        config = copy.copy(trackData)
        clientConfig = trackData["style"]
        del config["style"]

        if "match" in gffOpts:
            config["glyph"] = "JBrowse/View/FeatureGlyph/Segments"
            cmd += ["--type", gffOpts["match"]]

        cmd += ["--clientConfig", json.dumps(clientConfig)]

        if "trackType" in gffOpts:
            cmd += ["--trackType", gffOpts["trackType"]]
            if gffOpts["trackType"] == "BlastView/View/Track/CanvasFeatures":
                config["glyph"] = "JBrowse/View/FeatureGlyph/Segments"
        else:
            cmd += ["--trackType", "JBrowse/View/Track/CanvasFeatures"]

        if metadata:
            config.update({"metadata": metadata})
        cmd.extend(["--config", json.dumps(config)])

        self.subprocess_check_call(cmd)

    def add_rest(self, url, trackData):
        data = {
            "label": trackData["label"],
            "key": trackData["key"],
            "category": trackData["category"],
            "type": "JBrowse/View/Track/HTMLFeatures",
            "storeClass": "JBrowse/Store/SeqFeature/REST",
            "baseUrl": url,
            "query": {"organism": "tyrannosaurus"},
        }
        self._add_track_json(data)

    def process_annotations(self, track):
        category = track["category"].replace("__pd__date__pd__", TODAY)
        outputTrackConfig = {
            "style": {
                "label": track["style"].get("label", "description"),
                "className": track["style"].get("className", "feature"),
                "description": track["style"].get("description", ""),
            },
            "overridePlugins": track["style"].get("overridePlugins", False) == "True",
            "overrideDraggable": track["style"].get("overrideDraggable", False)
            == "True",
            "maxHeight": track["style"].get("maxHeight", "600"),
            "category": category,
        }

        for (
            i,
            (dataset_path, dataset_ext, track_human_label, extra_metadata),
        ) in enumerate(track["trackfiles"]):
            log.info("Processing %s / %s", category, track_human_label)
            outputTrackConfig["key"] = track_human_label
            # We add extra data to hash for the case of REST + SPARQL.
            try:
                rest_url = track["conf"]["options"]["url"]
            except KeyError:
                rest_url = ""

            # I chose to use track['category'] instead of 'category' here. This
            # is intentional. This way re-running the tool on a different date
            # will not generate different hashes and make comparison of outputs
            # much simpler.
            hashData = [dataset_path, track_human_label, track["category"], rest_url]
            hashData = "|".join(hashData).encode("utf-8")
            outputTrackConfig["label"] = hashlib.md5(hashData).hexdigest() + "_%s" % i

            # Colour parsing is complex due to different track types having
            # different colour options.
            colourOptions = self.cs.parse_colours(
                track["conf"]["options"], track["format"], gff3=dataset_path
            )
            # This used to be done with a dict.update() call, however that wiped out any previous style settings...
            for key in colourOptions:
                if key == "style":
                    for subkey in colourOptions["style"]:
                        outputTrackConfig["style"][subkey] = colourOptions["style"][
                            subkey
                        ]
                else:
                    outputTrackConfig[key] = colourOptions[key]

            # import sys; sys.exit()
            if dataset_ext in ("gff", "gff3", "bed"):
                self.add_features(
                    dataset_path,
                    dataset_ext,
                    outputTrackConfig,
                    track["conf"]["options"]["gff"],
                    metadata=extra_metadata,
                )
            elif dataset_ext == "bigwig":
                self.add_bigwig(
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["wiggle"],
                    metadata=extra_metadata,
                )
            elif dataset_ext == "bam":
                real_indexes = track["conf"]["options"]["pileup"]["bam_indices"][
                    "bam_index"
                ]
                if not isinstance(real_indexes, list):
                    # <bam_indices>
                    #  <bam_index>/path/to/a.bam.bai</bam_index>
                    # </bam_indices>
                    #
                    # The above will result in the 'bam_index' key containing a
                    # string. If there are two or more indices, the container
                    # becomes a list. Fun!
                    real_indexes = [real_indexes]

                self.add_bam(
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["pileup"],
                    bam_index=real_indexes[i],
                    metadata=extra_metadata,
                )
            elif dataset_ext == "blastxml":
                self.add_blastxml(
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["blast"],
                    metadata=extra_metadata,
                )
            elif dataset_ext == "vcf":
                self.add_vcf(dataset_path, outputTrackConfig, metadata=extra_metadata)
            elif dataset_ext == "rest":
                self.add_rest(
                    track["conf"]["options"]["url"],
                    outputTrackConfig,
                    metadata=extra_metadata,
                )
            else:
                log.warn("Do not know how to handle %s", dataset_ext)

            # Return non-human label for use in other fields
            yield outputTrackConfig["label"]

    def add_final_data(self, data):
        viz_data = {}
        if len(data["visibility"]["default_on"]) > 0:
            viz_data["defaultTracks"] = ",".join(data["visibility"]["default_on"])

        if len(data["visibility"]["always"]) > 0:
            viz_data["alwaysOnTracks"] = ",".join(data["visibility"]["always"])

        if len(data["visibility"]["force"]) > 0:
            viz_data["forceTracks"] = ",".join(data["visibility"]["force"])

        generalData = {}
        if data["general"]["aboutDescription"] is not None:
            generalData["aboutThisBrowser"] = {
                "description": data["general"]["aboutDescription"].strip()
            }

        generalData["view"] = {"trackPadding": data["general"]["trackPadding"]}
        generalData["shareLink"] = data["general"]["shareLink"] == "true"
        generalData["show_tracklist"] = data["general"]["show_tracklist"] == "true"
        generalData["show_nav"] = data["general"]["show_nav"] == "true"
        generalData["show_overview"] = data["general"]["show_overview"] == "true"
        generalData["show_menu"] = data["general"]["show_menu"] == "true"
        generalData["hideGenomeOptions"] = (
            data["general"]["hideGenomeOptions"] == "true"
        )
        generalData["plugins"] = data["plugins"]

        viz_data.update(generalData)
        self._add_json(viz_data)

        if "GCContent" in data["plugins_python"]:
            self._add_track_json(
                {
                    "storeClass": "JBrowse/Store/SeqFeature/SequenceChunks",
                    "type": "GCContent/View/Track/GCContentXY",
                    "label": "GCContentXY",
                    "urlTemplate": "seq/{refseq_dirpath}/{refseq}-",
                    "bicolor_pivot": 0.5
                    # TODO: Expose params for everyone.
                }
            )

        if "ComboTrackSelector" in data["plugins_python"]:
            with open(
                os.path.join(self.outdir, "data", "trackList.json"), "r"
            ) as handle:
                trackListJson = json.load(handle)
                trackListJson.update(
                    {
                        "trackSelector": {
                            "renameFacets": {
                                "tool_tool": "Tool ID",
                                "tool_tool_id": "Tool ID",
                                "tool_tool_version": "Tool Version",
                                "dataset_edam_format": "EDAM",
                                "dataset_size": "Size",
                                "history_display_name": "History Name",
                                "history_user_email": "Owner",
                                "metadata_dbkey": "Dbkey",
                            },
                            "displayColumns": [
                                "key",
                                "tool_tool",
                                "tool_tool_version",
                                "dataset_edam_format",
                                "dataset_size",
                                "history_display_name",
                                "history_user_email",
                                "metadata_dbkey",
                            ],
                            "type": "Faceted",
                            "title": ["Galaxy Metadata"],
                            "escapeHTMLInData": False,
                        },
                        "trackMetadata": {
                            "indexFacets": [
                                "category",
                                "key",
                                "tool_tool_id",
                                "tool_tool_version",
                                "dataset_edam_format",
                                "history_user_email",
                                "history_display_name",
                            ]
                        },
                    }
                )
                with open(
                    os.path.join(self.outdir, "data", "trackList2.json"), "w"
                ) as handle:
                    json.dump(trackListJson, handle)

    def clone_jbrowse(self, jbrowse_dir, destination):
        """Clone a JBrowse directory into a destination directory.
        """
        # JBrowse seems to have included some bad symlinks, cp ignores bad symlinks
        # unlike copytree
        cmd = ["cp", "-r", os.path.join(jbrowse_dir, "."), destination]
        log.debug(" ".join(cmd))
        subprocess.check_call(cmd)
        cmd = ["mkdir", "-p", os.path.join(destination, "data", "raw")]
        log.debug(" ".join(cmd))
        subprocess.check_call(cmd)

        # http://unix.stackexchange.com/a/38691/22785
        # JBrowse releases come with some broken symlinks
        cmd = ["find", destination, "-type", "l", "-xtype", "l"]
        log.debug(" ".join(cmd))
        symlinks = subprocess.check_output(cmd)
        for i in symlinks:
            try:
                os.unlink(i)
            except:
                pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("xml", type=argparse.FileType("r"), help="Track Configuration")

    parser.add_argument("--jbrowse", help="Folder containing a jbrowse release")
    parser.add_argument("--outdir", help="Output directory", default="out")
    parser.add_argument(
        "--standalone",
        help="Standalone mode includes a copy of JBrowse",
        action="store_true",
    )
    args = parser.parse_args()

    tree = ET.parse(args.xml.name)
    root = tree.getroot()

    jc = JbrowseConnector(
        jbrowse=args.jbrowse,
        outdir=args.outdir,
        genomes=[
            {
                "path": os.path.realpath(x.attrib["path"]),
                "meta": metadata_from_node(x.find("metadata")),
            }
            for x in root.findall("metadata/genomes/genome")
        ],
        standalone=args.standalone,
        gencode=root.find("metadata/gencode").text,
    )

    extra_data = {
        "visibility": {"default_on": [], "default_off": [], "force": [], "always": []},
        "general": {
            "defaultLocation": root.find("metadata/general/defaultLocation").text,
            "trackPadding": int(root.find("metadata/general/trackPadding").text),
            "shareLink": root.find("metadata/general/shareLink").text,
            "aboutDescription": root.find("metadata/general/aboutDescription").text,
            "show_tracklist": root.find("metadata/general/show_tracklist").text,
            "show_nav": root.find("metadata/general/show_nav").text,
            "show_overview": root.find("metadata/general/show_overview").text,
            "show_menu": root.find("metadata/general/show_menu").text,
            "hideGenomeOptions": root.find("metadata/general/hideGenomeOptions").text,
        },
        "plugins": [],
        "plugins_python": [],
    }

    plugins = root.find("plugins").attrib
    if plugins["GCContent"] == "True":
        extra_data["plugins_python"].append("GCContent")
        extra_data["plugins"].append(
            {
                "location": "https://cdn.rawgit.com/elsiklab/gccontent/5c8b0582ecebf9edf684c76af8075fb3d30ec3fa/",
                "name": "GCContent",
            }
        )

    if plugins["Bookmarks"] == "True":
        extra_data["plugins"].append(
            {
                "location": "https://cdn.rawgit.com/TAMU-CPT/bookmarks-jbrowse/5242694120274c86e1ccd5cb0e5e943e78f82393/",
                "name": "Bookmarks",
            }
        )

    if plugins["ComboTrackSelector"] == "True":
        extra_data["plugins_python"].append("ComboTrackSelector")
        extra_data["plugins"].append(
            {
                "location": "https://cdn.rawgit.com/TAMU-CPT/ComboTrackSelector/07f4a2d3434e7c8fd1c4cfa24c93b1e07c03029f/",
                "icon": "https://pbs.twimg.com/profile_images/580742252043100161/x64IwBFv_normal.png",
                "name": "ComboTrackSelector",
            }
        )

    if plugins["theme"] == "Minimalist":
        extra_data["plugins"].append(
            {
                "location": "https://cdn.rawgit.com/erasche/jbrowse-minimalist-theme/d698718442da306cf87f033c72ddb745f3077775/",
                "name": "MinimalistTheme",
            }
        )
    elif plugins["theme"] == "Dark":
        extra_data["plugins"].append(
            {
                "location": "https://cdn.rawgit.com/erasche/jbrowse-dark-theme/689eceb7e33bbc1b9b15518d45a5a79b2e5d0a26/",
                "name": "DarkTheme",
            }
        )

    # GALAXY_INFRASTRUCTURE_URL = root.find('metadata/galaxyUrl').text
    GALAXY_INFRASTRUCTURE_URL = "https://cpt.tamu.edu/galaxy/"
    for track in root.findall("tracks/track"):
        track_conf = {}
        track_conf["trackfiles"] = []

        for x in track.findall("files/trackFile"):
            metadata = metadata_from_node(x.find("metadata"))

            track_conf["trackfiles"].append(
                (
                    os.path.realpath(x.attrib["path"]),
                    x.attrib["ext"],
                    x.attrib["label"],
                    metadata,
                )
            )

        track_conf["category"] = track.attrib["cat"]
        track_conf["format"] = track.attrib["format"]
        try:
            # Only pertains to gff3 + blastxml. TODO?
            track_conf["style"] = {t.tag: t.text for t in track.find("options/style")}
        except TypeError as te:
            track_conf["style"] = {}
            pass
        track_conf["conf"] = etree_to_dict(track.find("options"))
        keys = jc.process_annotations(track_conf)

        for key in keys:
            extra_data["visibility"][
                track.attrib.get("visibility", "default_off")
            ].append(key)

    jc.add_final_data(extra_data)
