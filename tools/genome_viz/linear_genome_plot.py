#!/usr/bin/env python
from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator, GraphicRecord
from matplotlib import rc_context
import matplotlib
import matplotlib.pyplot as plt
from itertools import cycle
import re
import sys
import argparse

class CPTTranslator(BiopythonTranslator):
    """
    This is a customized translator from the dna_features_viewer module to fit Galaxy
    """

    global custom_feature_colors
    global box_status
    global label_fields
    global custom_name_colors
    global ignored_features_types
    global ignored_gene_labels
    global ignored_feature_labels

    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            if "product" in feature.qualifiers:
                color_specific = any(re.search(("(\\b"+str(item)+"\\b)"),feature.qualifiers["product"][0]) for item in custom_name_colors.keys()) or any(re.search((item),feature.qualifiers["product"][0]) for item in custom_name_colors.keys())
                if color_specific:
                    try:
                        return custom_name_colors[feature.qualifiers["product"][0]]
                    except KeyError:
                        for item in custom_name_colors.keys():
                            if item in feature.qualifiers["product"][0]:
                                custom_name_colors[feature.qualifiers["product"][0]] = custom_name_colors[item]
                                return custom_name_colors[feature.qualifiers["product"][0]]
                                #print(feature.qualifiers["product"][0])
                else:
                    try:
                        return custom_feature_colors[feature.type]
                    except KeyError:
                        return BiopythonTranslator.compute_feature_color(self, feature)
        else:
            if feature.type not in ignored_features_types:
                try:
                    return custom_feature_colors[feature.type]
                except KeyError:
                    return BiopythonTranslator.compute_feature_color(self, feature)

    def compute_feature_label(self, feature): # remove the chop_blocks
        self.label_fields = label_fields
        if feature.type == "CDS":
            if "product" in feature.qualifiers:
                if ignored_gene_labels: #  product name drop
                    verify_chops = any(re.search(("(\\b"+str(item)+"\\b)"),feature.qualifiers["product"][0]) for item in ignored_gene_labels) or any(re.search((item), feature.qualifiers["product"][0]) for item in ignored_gene_labels)
                    if verify_chops:
                        return None
                    else:
                        return BiopythonTranslator.compute_feature_label(self, feature)
                else:
                    return BiopythonTranslator.compute_feature_label(self, feature)
        elif feature.type in ignored_feature_labels:
            return None
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)

    def compute_filtered_features(self, features):
        return [
            feature for feature in features if feature.type not in ignored_features_types
        ]

    def compute_feature_legend_text(self, feature):
        return feature.type

    def compute_feature_box_color(self, feature):
        if feature.type == "CDS":
            return "white"

    def compute_feature_label_link_color(self, feature):
        return "black"

    def compute_feature_box_linewidth(self, feature):
        if box_status:
            return 0.5
        else:
            return 0

def parse_gbk(file):
    """ simple function to parse out the feature information AND products """

    record = SeqIO.read(file,"genbank")
    count = 0
    feature_types = {}
    product_names = []
    for feat in record.features:
        if feat.type not in feature_types:
            feature_types[feat.type] = 1
        else:
            feature_types[feat.type] += 1
        if "product" in feat.qualifiers:
            product_names.append(feat.qualifiers["product"][0])

    return feature_types, product_names, record

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Linear Genome Plot")
    #  Input and Parameters
    parser.add_argument("input_file",type=argparse.FileType("r"),help="genbank or gff3 file")
    parser.add_argument("--plot_width",type=int,default=20)
    #parser.add_argument("--plot_height",type=int,default=4)
    parser.add_argument("--title",type=str,default="genome plot") # NEED TO ADD TO XML
    parser.add_argument("--common_features_excluded", default="", help="common features to be excluded")
    parser.add_argument("--features_excluded",default="",help="features to be excluded from plot, separate by commas")
    parser.add_argument("--common_ignore_feature_labels", default="", help="common feature labels to be excluded")
    parser.add_argument("--ignored_feature_labels",default="",help="ignore labeling of specific features")
    parser.add_argument("--common_ignore_product_labels", default="", help="common product names to not label")
    parser.add_argument("--ignore_labeling",default="",help="labeling for specific products to ignore, separate by commas")
    parser.add_argument("--feature_label_order",default="locus_tag",help="label order, where the first choice is the first feature listed to pull name labels from") # NEED TO ADD TO XML
    parser.add_argument("--label_box",action="store_true",help="Use to have label box around feature labels")
    parser.add_argument("--label_algo",action="store_true",help="use dna features spacing algo for label placement (in or above feature)")
    #parser.add_argument("--level_offset",type=int,default=0,help="All features and annotations will be pushed up by the input amount. Useful for when plotting several sets of features successively on the same axis.") # Will exclude for now
    #parser.add_argument("--custom_region",action="store_true",help="cropped region for plot")
    parser.add_argument("--sz",type=int,help="beginning location for crop")
    parser.add_argument("--ez",type=int,help="end location for crop")
    parser.add_argument("--st",type=int,help="start site of translation")
    parser.add_argument("--et",type=int,help="end site of translation")
    parser.add_argument("--translation_on",action="store_true",help="plot the translation sub-axis")
    parser.add_argument("--feature_id",nargs="*",action="append",help="feature label to have custom color") # NEED TO ADD TO XML
    parser.add_argument("--feature_id_color",nargs="*",action="append",help="feature's accompanying color")
    parser.add_argument("--gene_id",nargs="*",action="append",help="gene/cds label to have custom color")
    parser.add_argument("--gene_id_color",nargs="*",action="append",help="gene/cds's accompanying color")
    parser.add_argument("--multiline",action="store_true",help="Plot multiline plot")
    parser.add_argument("--nucl_per_line",type=int,help="nucleotides per line of multiline")
    #  Output
    parser.add_argument("--file_stats",type=argparse.FileType("w"),default="out_stats.txt",help="output stat file")
    #parser.add_argument("--tmp_img",dest="tmp_img",type=argparse.FileType("wb"),default="out_tmp.svg")
    parser.add_argument("--out_img",dest="out_img",type=argparse.FileType("wb"),default="out_img.svg",help="svg genome plot")
    args = parser.parse_args()


    ##  Part I ; Parse and send output of features count and the list of product names
    feature_counts, products, genome = parse_gbk(args.input_file)
    with args.file_stats as f:
        f.writelines("---::: FILE BREAKDOWN :::---\n\n")
        f.writelines("------::: Feature Count :::------\n")
        for feature, count in feature_counts.items():
            f.writelines(f"Feature: {feature} ::::: Count: {count}\n")
        f.writelines("------::: Product Names :::------\n")
        if products != []:
            for product in products:
                f.writelines(f"Product Name: {product}\n")
        else:
            f.writelines("No Annotated Product Names Found")

    ##  Part II ; Prep Global Variables
    ##  Make K:V pairs for Feature Colors
    if args.label_box:
        box_status = True
    else:
        box_status = False

    if args.feature_id:
        feature_ids = [f for listed_obj in args.feature_id for f in listed_obj]
        feature_ids_colors = [f for listed_obj in args.feature_id_color for f in listed_obj]
        custom_feature_colors = dict(zip(feature_ids,feature_ids_colors))
    else:
        custom_feature_colors = {}

    ##  Make K:V pairs for Name Colors (as above)
    if args.gene_id:
        gene_ids = [g for listed_obj in args.gene_id for g in listed_obj]
        gene_ids_colors = [g for listed_obj in args.gene_id_color for g in listed_obj]
        custom_name_colors = dict(zip(gene_ids,gene_ids_colors))
    else:
        custom_name_colors = {}

    ##  Ignored Features
    #ignored_features_types = str.split(args.features_excluded,",")
    if args.common_features_excluded:
        ignored_features_types = str.split(args.common_features_excluded, ",")
        if args.features_excluded:
            ignored_features_types += str.split(args.features_excluded,",")
    elif args.features_excluded:
        ignored_features_types = str.split(args.features_excluded,",")
    else:
        ignored_features_types = False

    print(ignored_features_types)
    
    ## product labels
    if args.common_ignore_product_labels:
        ignored_gene_labels = str.split(args.common_ignore_product_labels,",")
        if args.ignore_labeling:
            ignored_gene_labels += str.split(args.ignore_labeling,",")
    elif args.ignore_labeling:
        ignored_gene_labels = str.split(args.ignore_labeling,",")
    else:
        ignored_gene_labels = False
    
    print(ignored_gene_labels)

    if args.feature_label_order != ['']:
        label_fields = str.split(args.feature_label_order,",")

    #if ignored_gene_labels == ['']:
    #    ignored_gene_labels = False

    ##  Ignored Labeling
    if args.common_ignore_feature_labels:
        ignored_feature_labels = str.split(args.common_ignore_feature_labels,",")
        if args.ignored_feature_labels:
            ignored_feature_labels += str.split(args.ignored_feature_labels,",")
    elif args.ignored_feature_labels:
        ignored_feature_labels = str.split(args.ignored_feature_labels,",")
    else:
        ignored_feature_labels = False
    
    print(ignored_feature_labels)
    ##  Print Statements for Debugging
    #print(custom_feature_colors)
    #print(custom_name_colors)
    #print(ignored_features_types)
    #print(ignored_gene_labels)
    #print(label_fields)

    ## Part III ; PLOT
    # Housekeeping
    rc_context({"font.family": ["monospace"],}) # courier-like
    matplotlib.use('Agg') # I think this has to be used...

    if args.label_algo:
        lab_algo = True
    else:
        lab_algo = False

    translator = CPTTranslator()
    graphic_record = translator.translate_record(genome)

    with open("tmp.svg", "wb") as img:
        img.truncate(0)
        img.close()

    if args.sz and not args.multiline: #  if user is wanting to look at a subset region of the genome
        zoom_start, zoom_end = args.sz, args.ez
        cropped = graphic_record.crop((zoom_start,zoom_end))
        ax, _ = cropped.plot(figure_width=args.plot_width, annotate_inline=lab_algo,figure_height=None)
        if args.translation_on:
            crop_seq = (args.st - 1, args.et)
            cropped.plot_translation(ax, location=crop_seq, fontdict={'size':8, 'weight':'bold'},y_offset=1)
        ax.set_title(args.title)
        # Galaxy specific shenanigans
        tmp_fig = "./tmp.svg"
        plt.savefig(tmp_fig)
        plt.close()
    elif args.multiline:
        if args.sz:
            zoom_start, zoom_end = args.sz, args.ez
        else:
            zoom_start, zoom_end = 1, graphic_record.sequence_length
        cropped = graphic_record.crop((zoom_start,zoom_end))
        ax, _ = cropped.plot_on_multiple_lines(
            figure_width=args.plot_width,
            annotate_inline=lab_algo,
            figure_height=None,
            nucl_per_line=args.nucl_per_line,
            plot_sequence=False
        )
        #ax.set_title(args.title)
        tmp_fig = "./tmp.svg"
        plt.savefig(tmp_fig)
        plt.close()
    else:
        ax, _ = graphic_record.plot(figure_width=args.plot_width, annotate_inline=lab_algo)
        ax.set_title(args.title)
        tmp_fig = "./tmp.svg"
        # Galaxy specific shenanigans
        plt.savefig(tmp_fig)
        plt.close()
    with open("tmp.svg", "rb") as img:
        for line in img:
            args.out_img.write(line)
