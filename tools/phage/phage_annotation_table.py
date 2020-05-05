#!/usr/bin/env python
# vim: set fileencoding=utf-8
import os
import argparse
from gff3 import genes, get_gff3_id, get_rbs_from, feature_test_true, feature_lambda, feature_test_type
from BCBio import GFF
from Bio import SeqIO
from jinja2 import Environment, FileSystemLoader
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="pat")

# Path to script, required because of Galaxy.
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
# Path to the HTML template for the report

def genes_all(feature_list, feature_type=["gene"], sort=False):
    """
    Simple filter to extract gene features from the feature set.
    """

    if not sort:
        for x in feature_lambda(
            feature_list, feature_test_type, {"types": feature_type}, subfeatures=True
        ):
            yield x
    else:
        data = list(genes_all(feature_list, feature_type, sort=False))
        data = sorted(data, key=lambda feature: feature.location.start)
        for x in data:
            yield x


def annotation_table_report(record, types, wanted_cols, gaf_data):
    getTypes = []
    for x in [y.strip() for y in types.split(",")]:
        getTypes.append(x)
    getTypes.append("gene")
    sorted_features = list(genes_all(record.features, getTypes, sort=True))
    if wanted_cols is None or len(wanted_cols.strip()) == 0:
        return [], []

    def rid(record, feature):
        """Organism ID
        """
        return record.id

    def id(record, feature):
        """ID
        """
        return feature.id

    def featureType(record, feature):
        """Type
        """
        return feature.type

    def name(record, feature):
        """Name
        """
        return feature.qualifiers.get("Name", ["None"])[0]

    def start(record, feature):
        """Boundary
        """
        return str(feature.location.start + 1)

    def end(record, feature):
        """Boundary
        """
        return str(feature.location.end)

    def location(record, feature):
        """Location
        """
        return str(feature.location.start + 1) + "..{0.end}".format(feature.location)

    def length(record, feature):
        """CDS Length (AA)
        """
        cdss = list(genes(feature.sub_features, feature_type="CDS", sort=True))
        return str(sum([len(cds) for cds in cdss]) / 3)

    def notes(record, feature):
        """User entered Notes"""
        return feature.qualifiers.get("Note", [])

    def date_created(record, feature):
        """Created"""
        return feature.qualifiers.get("date_creation", ["None"])[0]

    def date_last_modified(record, feature):
        """Last Modified"""
        return feature.qualifiers.get("date_last_modified", ["None"])[0]

    def description(record, feature):
        """Description"""
        return feature.qualifiers.get("description", ["None"])[0]

    def owner(record, feature):
        """Owner

        User who created the feature. In a 464 scenario this may be one of
        the TAs."""
        return feature.qualifiers.get("owner", ["None"])[0]

    def product(record, feature):
        """Product

        User entered product qualifier (collects "Product" and "product"
        entries)"""
        return feature.qualifiers.get(
            "product", feature.qualifiers.get("Product", ["None"])
        )[0]

    def note(record, feature):
        """Note

        User entered Note qualifier(s)"""
        return feature.qualifiers.get("Note", [])

    def strand(record, feature):
        """Strand
        """
        return "+" if feature.location.strand > 0 else "-"

    def sd_spacing(record, feature):
        """Shine-Dalgarno spacing
        """
        rbss = get_rbs_from(gene)
        if len(rbss) == 0:
            return "None"
        else:
            resp = []
            for rbs in rbss:
                cdss = list(genes(feature.sub_features, feature_type="CDS", sort=True))

                if rbs.location.strand > 0:
                    distance = min(
                        cdss, key=lambda x: x.location.start - rbs.location.end
                    )
                    distance_val = str(distance.location.start - rbs.location.end)
                    resp.append(distance_val)
                else:
                    distance = min(
                        cdss, key=lambda x: x.location.end - rbs.location.start
                    )
                    distance_val = str(rbs.location.start - distance.location.end)
                    resp.append(distance_val)

            if len(resp) == 1:
                return str(resp[0])
            return resp

    def sd_seq(record, feature):
        """Shine-Dalgarno sequence
        """
        rbss = get_rbs_from(gene)
        if len(rbss) == 0:
            return "None"
        else:
            resp = []
            for rbs in rbss:
                resp.append(str(rbs.extract(record).seq))
            if len(resp) == 1:
                return str(resp[0])
            else:
                return resp

    def start_codon(record, feature):
        """Start Codon
        """
        cdss = list(genes(feature.sub_features, feature_type="CDS", sort=True))
        data = [x for x in cdss]
        if len(data) == 1:
            return str(data[0].extract(record).seq[0:3])
        else:
            return [
                "{0} ({1.location.start}..{1.location.end}:{1.location.strand})".format(
                    x.extract(record).seq[0:3], x
                )
                for x in data
            ]

    def stop_codon(record, feature):
        """Stop Codon
        """
        return str(feature.extract(record).seq[-3:])

    def dbxrefs(record, feature):
        """DBxrefs
        """
        return feature.qualifiers.get("Dbxref", [])

    def upstream_feature(record, feature):
        """Next gene upstream"""
        if feature.strand > 0:
            upstream_features = [
                x for x in sorted_features if (x.location.start < feature.location.start and x.type == "gene" and x.strand == feature.strand)
            ]
            if len(upstream_features) > 0:
                return upstream_features[-1]
            else:
                return None
        else:
            upstream_features = [
                x for x in sorted_features if (x.location.end > feature.location.end and x.type == "gene" and x.strand == feature.strand)
            ]

            if len(upstream_features) > 0:
                return upstream_features[0]
            else:
                return None

    def up_feat(record, feature):
        """Next gene upstream"""
        up = upstream_feature(record, feature)
        if up:
            return str(up.id)
        return "None"

    def ig_dist(record, feature):
        """Distance to next upstream gene on same strand"""
        up = upstream_feature(record, feature)
        if up:
            dist = None
            if feature.strand > 0:
                dist = feature.location.start - up.location.end
            else:
                dist = up.location.start - feature.location.end
            return str(dist)
        else:
            return "None"

    def _main_gaf_func(record, feature, gaf_data, attr):
        if feature.id in gaf_data:
            return [x[attr] for x in gaf_data[feature.id]]
        return []

    def gaf_annotation_extension(record, feature, gaf_data):
        """GAF Annotation Extension

        Contains cross references to other ontologies that can be used
        to qualify or enhance the annotation. The cross-reference is
        prefaced by an appropriate GO relationship; references to
        multiple ontologies can be entered. For example, if a gene
        product is localized to the mitochondria of lymphocytes, the GO
        ID (column 5) would be mitochondrion ; GO:0005439, and the
        annotation extension column would contain a cross-reference to
        the term lymphocyte from the Cell Type Ontology.
        """
        return _main_gaf_func(record, feature, gaf_data, "annotation_extension")

    def gaf_aspect(record, feature, gaf_data):
        """GAF Aspect code

        E.g. P (biological process), F (molecular function) or C (cellular component)
        """
        return _main_gaf_func(record, feature, gaf_data, "aspect")

    def gaf_assigned_by(record, feature, gaf_data):
        """GAF Creating Organisation
        """
        return _main_gaf_func(record, feature, gaf_data, "assigned_by")

    def gaf_date(record, feature, gaf_data):
        """GAF Creation Date
        """
        return _main_gaf_func(record, feature, gaf_data, "date")

    def gaf_db(record, feature, gaf_data):
        """GAF DB
        """
        return _main_gaf_func(record, feature, gaf_data, "db")

    def gaf_db_reference(record, feature, gaf_data):
        """GAF DB Reference
        """
        return _main_gaf_func(record, feature, gaf_data, "db_reference")

    def gaf_evidence_code(record, feature, gaf_data):
        """GAF Evidence Code
        """
        return _main_gaf_func(record, feature, gaf_data, "evidence_code")

    def gaf_go_id(record, feature, gaf_data):
        """GAF GO ID
        """
        return _main_gaf_func(record, feature, gaf_data, "go_id")

    def gaf_go_term(record, feature, gaf_data):
        """GAF GO Term
        """
        return _main_gaf_func(record, feature, gaf_data, "go_term")

    def gaf_id(record, feature, gaf_data):
        """GAF ID
        """
        return _main_gaf_func(record, feature, gaf_data, "id")

    def gaf_notes(record, feature, gaf_data):
        """GAF Notes
        """
        return _main_gaf_func(record, feature, gaf_data, "notes")

    def gaf_owner(record, feature, gaf_data):
        """GAF Creator
        """
        return _main_gaf_func(record, feature, gaf_data, "owner")

    def gaf_with_or_from(record, feature, gaf_data):
        """GAF With/From
        """
        return _main_gaf_func(record, feature, gaf_data, "with_or_from")

    cols = []
    data = []
    funcs = []
    lcl = locals()
    for x in [y.strip().lower() for y in wanted_cols.split(",")]:
        if not x:
            continue
        if x == "type":
          x = "featureType"
        if x in lcl:
            funcs.append(lcl[x])
            # Keep track of docs
            func_doc = lcl[x].__doc__.strip().split("\n\n")
            # If there's a double newline, assume following text is the
            # "help" and the first part is the "name". Generate empty help
            # if not provided
            if len(func_doc) == 1:
                func_doc += [""]
            cols.append(func_doc)
        elif "__" in x:
            chosen_funcs = [lcl[y] for y in x.split("__")]
            func_doc = [
                " of ".join(
                    [y.__doc__.strip().split("\n\n")[0] for y in chosen_funcs[::-1]]
                )
            ]
            cols.append(func_doc)
            funcs.append(chosen_funcs)

    for gene in genes_all(record.features, getTypes, sort=True):
        row = []
        for func in funcs:
            if isinstance(func, list):
                # If we have a list of functions, repeatedly apply them
                value = gene
                for f in func:
                    if value is None:
                        value = "None"
                        break

                    value = f(record, value)
            else:
                # Otherwise just apply the lone function
                if func.func_name.startswith("gaf_"):
                    value = func(record, gene, gaf_data)
                else:
                    value = func(record, gene)

            if isinstance(value, list):
                collapsed_value = ", ".join(value)
                value = [str(collapsed_value).decode("utf-8")]
            else:
                value = str(value).decode("utf-8")

            row.append(value)
        # print row
        data.append(row)
    return data, cols


def parseGafData(file):
    cols = []
    data = {}
    # '10d04a01-5ed8-49c8-b724-d6aa4df5a98d': {
    # 'annotation_extension': '',
    # 'aspect': '',
    # 'assigned_by': 'CPT',
    # 'date': '2017-05-04T16:25:22.161916Z',
    # 'db': 'UniProtKB',
    # 'db_reference': 'GO_REF:0000100',
    # 'evidence_code': 'ISA',
    # 'gene': '0d307196-833d-46e8-90e9-d80f7a041d88',
    # 'go_id': 'GO:0039660',
    # 'go_term': 'structural constituent of virion',
    # 'id': '10d04a01-5ed8-49c8-b724-d6aa4df5a98d',
    # 'notes': 'hit was putative minor structural protein',
    # 'owner': 'amarc1@tamu.edu',
    # 'with_or_from': 'UNIREF90:B2ZYZ7'
    # },
    for row in file:
        if row.startswith("#"):
            # Header
            cols = (
                row.strip().replace("# ", "").replace("GO Term", "go_term").split("\t")
            )
        else:
            line = row.strip().split("\t")
            tmp = dict(zip(cols, line))
            if tmp["gene"] not in data:
                data[tmp["gene"]] = []

            data[tmp["gene"]].append(tmp)
    return data


def evaluate_and_report(
    annotations,
    genome,
    types="gene",
    reportTemplateName="phage_annotation_validator.html",
    annotationTableCols="",
    gafData=None,
):
    """
    Generate our HTML evaluation of the genome
    """
    # Get features from GFF file
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    # Get the first GFF3 record
    # TODO: support multiple GFF3 files.
    at_table_data = []
    gaf = {}
    if gafData:
        gaf = parseGafData(gafData)

    for record in GFF.parse(annotations, base_dict=seq_dict):
        if reportTemplateName.endswith(".html"):
            record.id = record.id.replace(".", "-")
        log.info("Producing an annotation table for %s" % record.id)
        annotation_table_data, annotation_table_col_names = annotation_table_report(
            record, types, annotationTableCols, gaf
        )
        at_table_data.append((record, annotation_table_data))
        # break

    # This is data that will go into our HTML template
    kwargs = {
        "annotation_table_data": at_table_data,
        "annotation_table_col_names": annotation_table_col_names,
    }

    env = Environment(
        loader=FileSystemLoader(SCRIPT_PATH), trim_blocks=True, lstrip_blocks=True
    )
    if reportTemplateName.endswith(".html"):
        env.filters["nice_id"] = str(get_gff3_id).replace(".", "-")
    else:
        env.filters["nice_id"] = get_gff3_id

    def join(listy):
        return "\n".join(listy)

    env.filters.update({"join": join})
    tpl = env.get_template(reportTemplateName)
    return tpl.render(**kwargs).encode("utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="rebase gff3 features against parent locations", epilog=""
    )
    parser.add_argument(
        "annotations", type=argparse.FileType("r"), help="Parent GFF3 annotations"
    )
    parser.add_argument("genome", type=argparse.FileType("r"), help="Genome Sequence")

    parser.add_argument(
        "--types",
        help="Select extra types to display in output (Will always include gene)",
    )

    parser.add_argument(
        "--reportTemplateName",
        help="Report template file name",
        default="phageqc_report_full.html",
    )
    parser.add_argument(
        "--annotationTableCols",
        help="Select columns to report in the annotation table output format",
    )
    parser.add_argument(
        "--gafData", help="CPT GAF-like table", type=argparse.FileType("r")
    )

    args = parser.parse_args()

    print(evaluate_and_report(**vars(args)))
