from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator, GraphicRecord
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
    global label_fields
    global custom_name_colors
    global ignored_features_types
    global ignored_gene_labels

    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            if "product" in feature.qualifiers:
                color_specific = any(re.search(("(\\b"+str(item)+"\\b)"),feature.qualifiers["product"][0]) for item in custom_name_colors.keys())
                if color_specific:
                    print(custom_name_colors[feature.qualifiers["product"][0]])
                    return custom_name_colors[feature.qualifiers["product"][0]]
                else:
                    return custom_feature_colors[feature.type]
        else:
            if feature.type not in ignored_features_types:
                try:
                    return custom_feature_colors[feature.type]
                except KeyError:
                    sys.exit("ERROR: Features included for plotting do not match custom color schemea")

    def compute_feature_label(self, feature): # remove the chop_blocks
        self.label_fields = label_fields
        if feature.type == "CDS":
            if "product" in feature.qualifiers:
                verify_chops = any(re.search(("(\\b"+str(item)+"\\b)"),feature.qualifiers["product"][0]) for item in ignored_gene_labels)
                if verify_chops:
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
    
    def compute_featurebox_linewidth(self, feature):
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
    parser.add_argument("--title",type=str,default="genome plot") # NEED TO ADD TO XML
    parser.add_argument("--features_excluded",default="",help="features to be excluded from plot, separate by commas")
    parser.add_argument("--ignore_labeling",default="",help="labeling for specific genes to ignore, separate by commas")
    parser.add_argument("--feature_label_order",default="CDS",help="label order, where the first choice is the first feature listed to pull name labels from") # NEED TO ADD TO XML
    parser.add_argument("--label_above",action="store_true",help="force all labels above gene")
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
    #  Output
    parser.add_argument("--file_stats",type=argparse.FileType("w"),default="out_stats.txt",help="output stat file")
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
    ignored_features_types = str.split(args.features_excluded,",")
    ignored_gene_labels = str.split(args.ignore_labeling,",")
    label_fields = str.split(args.feature_label_order,",")

    ##  Print Statements for Debugging
    print(custom_feature_colors)
    print(custom_name_colors)
    print(ignored_features_types)
    print(ignored_gene_labels)

    ## Part III ; PLOT
    rc_context({"font.family": ["monospace"],})
        above = True
    else:
        above = False

    translator = CPTTranslator()
    graphic_record = translator.translate_record(genome)

    if args.sz: #  if user is wanting to look at a subset region of the genome
        print("-- crop mode --")
        zoom_start, zoom_end = args.sz, args.ez
        cropped = graphic_record.crop((zoom_start,zoom_end))
        ax, _ = cropped.plot(figure_width=args.plot_width, elevate_outline_annotations=True)
        if args.translation_on:
            crop_seq = (args.st - 1, args.et)
            cropped.plot_translation(ax, location=crop_seq, fontdict={'size':8, 'weight':'bold'},y_offset=1)
        ax.set_title(args.title)
        tmp_fig = "./tmp.svg"
        ax.figure.savefig(tmp_fig)
        plt.close()
        with open("tmp.svg", "rb") as img:
            for line in img:
                args.out_img.write(line)
        # do the chop cut
    else:
        ax, _ = graphic_record.plot(figure_width=args.plot_width, annotate_inline=above)
        ax.set_title(args.title)
        plt.close()
        tmp_fig = "./tmp.svg"
        ax.figure.savefig(tmp_fig)
        plt.close()
        with open("tmp.svg", "rb") as img:
            for line in img:
                args.out_img.write(line)
