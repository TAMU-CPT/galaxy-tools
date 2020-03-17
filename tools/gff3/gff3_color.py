#!/usr/bin/env python
import re
import sys
import logging
import argparse
from BCBio import GFF
from gff3 import feature_lambda, wa_unified_product_name

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


color_scheme = {
    "defense": {
        "members": [
            "rII",
            "rIIA",
            "rIIB",
            "rex",
            "rexA",
            "rexB",
            "ocr",
            "dar",
            "darA",
            "darB",
        ],
        "color": "#C8FFC8",
        "title": "Defense",
    },
    "trna": {"members": ["tRNA"], "color": "#71BC78", "title": "tRNAs"},
    "regulation": {
        "title": "Regulation",
        "color": "#FFA500",
        "members": [
            "Translational repressor",
            "RegA",
            "RegB",
            "regulatory",
            "regulator",
            "transcriptional repressor",
            "anti-repressor",
            "rna",
            "rna polymerase",
            "Sigma Factor",
        ],
    },
    "biosynthesis": {
        "title": "Biosynthesis",
        "color": "#0000FF",
        "members": [
            "Deoxynucleoside",
            "Deoxyribonucleotidase",
            "Deoxyuridine",
            "Ribonucleoside-diphoshate reductase",
            "Serine kinase",
            "threonine kinase",
            "cytidine deaminase",
            "dUTPase",
            "dUTPase",
            "deoxynucleotide",
            "dihydrofolate reductase",
            "glutaredoxin",
            "guanylate kinase",
            "reductase",
            "ribonucleotide reductase",
            "thioredoxin",
            "thymidylate",
        ],
    },
    "dna_rep_recomb": {
        "title": "DNA Replication/Recombination",
        "color": "#FFFF00",
        "members": [
            "Clamp",
            "DNA binding protein",
            # Should we be testing versions with .replace(' ', '').replace('-', '')?
            "DNA-binding protein",
            "DNA end Protector",
            "DNA ligase",
            "DexA",
            "DnaA",
            "DnaB",
            "DnaQ",
            "Helicase",
            "RNA ligase",
            "RNaseH",
            "RecA",
            "RecF",
            "Recombination",
            "RuvC",
            "UvsW",
            "UvsY",
            "helicase",
            "holliday junction",
            "phosphoesterase",
            "primase",
            "recombinase",
            "recombination",
            "repair",
            "single strand annealing",
            "topoisomerase",
            "whisker",
            "sliding",
            "methylase",
            "methyltransferase",
            "mom",
            "glucosyl\\s*transferase",
            "glycosyl\\s*transferase",
            "integrase",
        ],
        "custom": {
            "hnh": {"isnot": ["HNH", "homing endonuclease"], "is": ["nuclease"]},
            "polymerase": {
                "is": ["polymerase"],
                "isnot": ["rna polymerase", "polymerase sigma factor"],
            },
        },
    },
    "novel": {"members": ["Novel"], "color": "#AAAAAA", "title": "Novel"},
    "dna_pack": {
        "members": ["terminase"],
        "color": "#00FFFF",
        "title": "DNA Packaging",
    },
    "terminator": {
        "color": "#00FF00",
        "title": "terminator",
        "members": ["terminator"],
    },
    "lysis": {
        "custom": {
            "lysozyme": {
                "isnot": ["lysozyme baseplate", "tail lysozyme"],
                "is": ["lysozyme"],
            }
        },
        "members": [
            "antiholin",
            "holin",
            "endolysin",
            "spanin",
            "peptidoglycan",
            "amidase",
            "transglycosylase",
            "carboxypeptidase",
        ],
        "color": "#FF00FF",
        "title": "Lysis",
    },
    "morpho": {
        "custom": {
            "tail": {"isnot": ["tail lysozyme"], "is": ["tail"]},
            "baseplate": {"is": ["baseplate"], "isnot": ["lysozyme baseplate"]},
        },
        "members": [
            "tail\\s*spike",
            "fiber",
            "neck",
            "sheath",
            "tube",
            "pectin",
            "prohead",
            "scaffold",
            "capsid",
            "head",
            "head-to-tail joining",
            "pre-neck",
            "Tape",
            "tailspike",
            "structural",
            "morphogenesis",
            "assembly",
            "chaperone",
            "joining",
            "decoration",
            "protease",
            "frameshift",
            "portal",
        ],
        "color": "#87CEFA",
        "title": "Morphogenesis",
    },
    "hnh": {
        "title": "HNH/Homing/GIY-YIG",
        "color": "#C89664",
        "members": ["HNH", "homing endonuclease", "GIY-YIG"],
    },
}


class ColorScheme(object):
    def __init__(self):
        self.standard_regex = {}
        self.custom_regex = {}
        for key in color_scheme:
            for member in color_scheme[key]["members"]:
                regex = re.compile(member, re.IGNORECASE)
                color = color_scheme[key]["color"]
                self.standard_regex[member] = {
                    "str": member,
                    "regex": regex,
                    "color": color,
                }

            if "custom" in color_scheme[key]:
                for custom_key in color_scheme[key]["custom"]:
                    cu = color_scheme[key]["custom"][custom_key]
                    self.custom_regex[key + "/" + custom_key] = {
                        "color": color_scheme[key]["color"],
                        "is": [re.compile("\b" + x + "\b") for x in cu.get("is", [])],
                        "isnot": [
                            re.compile("\b" + x + "\b") for x in cu.get("isnot", [])
                        ],
                    }

    def get_color(self, product_list):

        for product in product_list:
            if product is None:
                continue

            matched = None
            for regex in self.standard_regex:

                if re.search(self.standard_regex[regex]["regex"], product):
                    matched = self.standard_regex[regex]["color"]

            for regex in self.custom_regex:
                care = False

                for re_is in self.custom_regex[regex]["is"]:
                    if re.search(re_is, product):
                        care = True

                for re_isnot in self.custom_regex[regex]["isnot"]:
                    if re.search(re_isnot, product):
                        care = False

                if care:
                    is_ok = True
                    ok_to_overwrite = True

                    # COPYPASTA from original perl:
                    # Fixes a strange bug. If we have a match, and then
                    # we have a subpart of that match which hits a custom
                    # element, we want to make sure that it'll ONLY overwrite
                    # if we didn't specifically exclude items like the one we
                    # hit.
                    for re_is in self.custom_regex[regex]["is"]:
                        if not re.search(re_is, product):
                            is_ok = False

                    for re_isnot in self.custom_regex[regex]["isnot"]:
                        if re.search(re_isnot, product):
                            ok_to_overwrite = False
                            is_ok = False

                    if is_ok and ok_to_overwrite:
                        matched = self.custom_regex[regex]["color"]

            if matched is not None:
                log.info("%s -> %s", product, matched)
                return matched


def apply_color(feature, **kwargs):
    product = [wa_unified_product_name(feature)]
    color = kwargs["cs"].get_color(product)
    if color is not None:
        feature.qualifiers["color"] = color
    return True


def gff_filter(gff3):
    cs = ColorScheme()

    for rec in GFF.parse(gff3):
        rec.features = feature_lambda(
            rec.features, apply_color, {"cs": cs}, subfeatures=False
        )
        rec.annotations = {}
        GFF.write([rec], sys.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="add color qualifiers based on product qualifiers"
    )
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 annotations")
    args = parser.parse_args()
    gff_filter(**vars(args))
