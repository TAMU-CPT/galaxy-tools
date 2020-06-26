from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator, GraphicRecord
import matplotlib.pyplot as plt
from itertools import cycle
import numpy as np
import argparse

class CPTTranslator(BiopythonTranslator):
    """ 
    This is a customized translator from the dna_features_viewer module to fit Galaxy
    """

    global ignored_features_types
    global ignored_gene_labels
    global label_fields
    global custom_name_colors
    global custom_feature_colors

    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            if "product" in feature.qualifiers:
                color_specific = any(re.search(("(\\b"+str(item)+"\\b)"),feature.qualifiers["product"][0]) for item in custom_name_colors.keys())
                if color_specific:
                    print(custom_name_colors[feature.qualifiers["product"][0]])
                    return custom_name_colors[feature.qualifiers["product"][0]]
                else:
                    return custom_feature_colors[feature.type]

    def compute_feature_label(self, feature): # remove the chop_blocks
        self.label_fields = label_fields
        if feature.type == "CDS":
            if "product" in feature.qualifiers:
                verify_chops = any(re.search(("(\\b"+str(item)+"\\b)"),feature.qualifiers["product"][0]) for item in chop_block)
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
    
    return feature_types, product_names

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Linear Genome Plot")
    #  Input and Parameters
    parser.add_argument("input_file",type=argparse.FileType("r"),help="genbank or gff3 file")
    parser.add_argument("--plot_width",type=int,default=20)
    parser.add_argument("--features_excluded",default="",help="features to be excluded from plot, separate by commas")
    parser.add_argument("--ignore_labeling",default="",help="labeling for specific genes to ignore, separate by commas")
    parser.add_argument("--label_above",action="store_true",help="force all labels above gene")
    #parser.add_argument("--custom_region",action="store_true",help="cropped region for plot")
    parser.add_argument("--sz",type=int,help="beginning location for crop")
    parser.add_argument("--ez",type=int,help="end location for crop")
    parser.add_argument("--st",type=int,help="start site of translation")
    parser.add_argument("--et",type=int,help="end site of translation")
    parser.add_argument("--feature_id",nargs="*",action="append",help="feature label to have custom color")
    parser.add_argument("--feature_id_color",nargs="*",action="append",help="feature's accompanying color")
    parser.add_argument("--gene_id",nargs="*",action="append",help="gene/cds label to have custom color")
    parser.add_argument("--gene_id_color",nargs="*",action="append",help="gene/cds's accompanying color")
    #  Output
    parser.add_argument("--file_stats",type=argparse.FileType("w"),default="out_stats.txt",help="output stat file")
    parser.add_argument("--out_img",type=argparse.FileType("w"),default="out_img.svg",help="svg genome plot")
    args = parser.parse_args()


    ##  Part I ; Parse and send output of features count and the list of product names
    feature_counts, products = parse_gbk(args.input_file)
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

    ##  Print Statements for Debugging
    print(custom_feature_colors)
    print(custom_name_colors)
    print(ignored_features_types)
    print(ignored_gene_labels)

    ## Part III ; PLOT

    #ignored_features_types
    #chop_block
    #label_fields
    #custom_name_colors
    #custom_feature_colors
